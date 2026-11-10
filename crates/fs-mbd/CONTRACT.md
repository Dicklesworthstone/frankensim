## The native rigid body as a composition-kernel Component

`component::RigidBody` wraps `RigidBodyIntegrator::step` behind
`fs_ports::component::Component` (contract v0), so that the kernel in
`fs-system` can join a free rigid body to other components without code
written for the combination. It adds no integrator and changes none.

Two ports, both receiving effort, three components each. `force` takes the
world-frame force on the centre of mass and returns the mean centre-of-mass
velocity `(x1 - x0) / dt`. `torque` takes the body-frame torque and returns the
mean body angular velocity `(omega0 + omega1) / 2`. These are the conjugate
increments the midpoint step already pairs with a held force and torque, so
the work the kernel books from the wire values equals `StepWorkDiagnostics::
applied_work_j` to roundoff. Stored energy is `DynamicsDiagnostics::
mechanical_energy`. The receipt's stored change is `mechanical_energy_change_j`
and its numerical loss is minus `work_energy_residual_j`; dissipation and
entropy production are zero. `RigidBodyRefusal` carries the native
`DynamicsError` unchanged, or `Cancelled` when the `Cx` was already requested.

The receipt's energy tolerance is 64 machine epsilons of three things added
together: the products of the norms of force and mean velocity and of torque
and mean angular velocity (the angular velocity comes out of a frame rotation,
so its components carry roundoff of the whole vector); the kinetic energies
and the sizes of the three terms of the gravity potential `-m g . x` at both
ends of the step (in free fall kinetic and potential energy cancel, and on a
level surface of an oblique field the three terms cancel each other, so the
measured total is only as exact as its terms); and the magnitudes of the
changes and of the residual. A tolerance scaled only by
the energy exchanged in the step refused free fall from the origin after 298
steps; that defect was found by an independent review on 2026-11-10 and is
pinned by the tests below. A second review the same day showed that the
magnitude of the potential's total was still too small under oblique gravity
far from the origin (the network refused at step 1, the single body failed
278 of 2,000 steps); both components now use the per-axis terms.

G0, in `crates/fs-system/tests/rigid_body_slice.rs` (run with
`cargo test -p fs-system --features slices`):
- a tumbling asymmetric body under force, torque and gravity passes
  `component::conformance::check_component`, and the work booked from the wire
  values equals the integrator's `applied_work_j` to 1e-13 of the paired
  magnitudes;
- 400 generated states pass the conformance check: heights to 500 m, fall
  speeds to 100 m/s, a principal frame rotated against the body frame, free
  fall, a torque across the spin, general wrenches, two step sizes;
- a body released at the origin under gravity follows the parabola for 5,000
  steps to 1e-10 with its stored energy unchanged to 1e-10 of the kinetic
  energy;
- the same body with zero torque, bonded by the kernel to a three-axis spring
  defined in the test, follows the monolithic midpoint rule for translation
  over 4,000 steps, keeps the size of its body angular momentum, and the
  window ledger closes.

No claim: no contact, constraint or collision; no analytic tangent of the port
map (a kernel loop differences it); the split variants of the integrator are
not wrapped; fs-mbd now depends on fs-ports for this module only.

## A damped rigid-body network as a composition-kernel Component

`component::network::ViscousNetwork` wraps
`step_viscous_network_with_loads` behind the same contract: any number of
bodies and passive viscous connectors, advanced by the owner's splitting (half
a connector source at fixed geometry, the native rigid motion of every body,
half a source). The wrapper integrates nothing.

Ports, in this order: `force<i>` per body (receives the world force on the
centre of mass, three components; returns the mean centre-of-mass velocity
`(x1 - x0) / dt`, the displacement the owner pairs with a held force, since
the sources act at fixed geometry) and `heat<e>` per connector (a heat port
that receives the temperature at which the connector's heat is released and
returns minus its mean dissipation rate, `-D_e / dt`). The state is the bank
of `RigidBodyState`. Stored energy is kinetic plus uniform-gravity potential.
The receipt's stored change, dissipation and connector heats are the owner's
numbers; the owner's work-heat residual is booked as numerical loss; entropy
is `sum(D_e / T_e)`. `ViscousNetworkRefusal` carries the owner's error
unchanged, or names a state of the wrong length, a temperature that is not
positive, or cancellation. `ViscousNetworkBuildError` refuses no bodies and a
connector that does not join two different bodies of the network.

G0, in `crates/fs-system/tests/viscous_network_slice.rs` (run with
`cargo test -p fs-system --features slices`): conformance under forces and
gravity at three step sizes with the kernel's booked energy equal to the
owner's work minus heat; a propose returns the owner's states bit for bit;
the refusals above; the damping-heat pairing as a declared graph with heat
routed by fraction to two thermal bodies and a named sink (one body shared by
both connectors), bit-identical to the owner marched alone over 400 steps;
held forces under gravity booked as boundary work. Three planted defects in
this module were each caught (half the heat leaving, velocities on the wrong
body, entropy at the wrong temperature).

No claim: no torque port, because the owner reports couple work as one total
and not per body, so a couple cannot be attached to a body through the kernel
yet; no analytic tangent; one propose allocates the load and state banks; the
owner's accuracy statements are unchanged by the wrapping; the thermal side
in the test is the kernel's ideal `HeatCapacity`, not the fs-thermal
deposition timeline that `fs-couple::rigid_damping_heat` drives.

## Affine friction stopping and reversal on finite bodies

The original rate-driven Hertz adapter retains the law's `exact_invariant`
flag and first stopping interval through every native trial and accepted body
publication. Admitted affine sliding and spinning paths can now stop, transmit
a supported static load, or reverse under an opposing load outside the original
static surface. The source curve, complete body mobility, force/couple transport
and held geometry remain the existing owners; no velocity threshold captures a
moving state.

Native source work includes common-body motion during a relative-contact hold.
Friction heat stops when relative motion stops, while the applied loads can
continue to accelerate the pair. Closed-form event receipts are distinguished
from embedded adaptive estimates; delegated curved Coulomb motion still needs
its independent timestep refinement. Focused analytic finite-body tests cover
translation/spin braking, reversal, static continuation, split replay, original
momenta and actual supplied work. Failed effort and cancellation gates retain
both bodies and the prior mode. General curved and non-affine stopping events
remain subject to the original law's explicit admission/refusal boundary.

## Driven original local-speed Hertz friction on finite bodies

`tangential_friction::hertz::rate_dependent::driven` advances the original
pressure-resolved friction curve under held world forces and body-frame torques.
The captured full contact mobility and native force/couple transport retain both
finite bodies. Every temperature/coefficient trial starts from the same original
momenta, external loads and explicit resting/sliding mode; load/radius changes
refuse. The constitutive initial-value trajectory follows the original weakening
branch without an algebraic terminal-root selection or centre-speed coefficient.

The private `DrivenRigidPair` is shared with the existing driven Coulomb owner.
It prepares external acceleration, applies one summed impulse per body, retains
common-body source work and independently checks actual momentum, kinetic change,
contact rates and work. Its closed-normal path gate consumes the complete rate
law's finite-speed/static wrench bound. Source work and physical friction heat
remain separate; unresolved driven zero-rate events publish neither bodies nor
mode. The original native report carries all constitutive numerical estimates
and their availability flag. Delegated Coulomb zeros do not assert a measured
zero time-discretization error.

Poses, contact arms, inertia, load, radius and external drive remain held within
a phase. Pressure resolution and adaptive time errors require independent
refinement; general curved or non-affine driven arrest, changing normal geometry
and source-regime inference remain outside this phase owner. Focused native regressions compare analytic
affine translation/spin/release and true full-body source work, an independent
mixed weakening body-coordinate ODE, world-frame transformation, reusable source
trials, geometry refusal, atomic mode publication and spectator preservation.

## Driven original-body Hertz shear and spin

`tangential_friction::hertz::driven` advances held external world forces and
body-frame torques together with the original joint Hertz sliding/spin law.
Original masses, full inertia, orientations and contact arms supply the complete
three-coordinate mobility. The explicit resting/sliding mode survives accepted
continuation. Static admission, breakaway and arrest belong to the original
pressure-resolved surface; no separate force/spin caps or velocity clamp replace it.

Each original body receives one summed external/contact impulse. Actual native
momenta and compensated kinetic change independently check the law endpoint and
work. External work retains common-body motion: the analytically integrated free
source work is corrected by `g dot (travel - free_travel)`, where `W g = a_drive`.
It is never inferred from heat, numerical loss or rounded free endpoint storage.
Static force transmission can accelerate both bodies with nonzero drive work and
zero frictional heat. The full balance is `Delta K + Q_friction - W_external`.
Bodies and retained mode publish together only after every original gate passes.

Poses, arms, inertia, normal load and radius remain held. A conservative full-path
normal-response bound refuses unsupported opening. Invariant directions retain
exact constant-acceleration events; curved motion and curved capture retain the
law's explicit first-order timestep approximation. Finite pose/normal evolution,
changing source drive and general driven rate-weakening laws remain separate.
Focused native tests cover analytic common-motion work, motor spin and split
continuation, an independent full-body power/trajectory ODE, world/body-frame
covariance, unrepresentable source work, atomic refusal and unaffected bodies.

## Native recoil from local-speed Hertz friction

`tangential_friction::hertz::rate_dependent` advances the original
`HertzFrictionPatch` against the captured finite body's complete contact mobility.
The actual translating/spinning pressure field supplies each local slip speed;
no centre-speed coefficient replaces it. The law's finite speed/coefficient and
whole-wrench bounds remain visible in the native report.

The original exact/curved pair owners share geometry preparation, one summed
point/couple impulse update and independent per-body momentum, endpoint and
energy audits. Preparation depends on geometry and inertia without constructing
a substitute Coulomb law. Repeated curve trials start from the same momenta and
retain the original normal load and radius. The whole-path normal-response gate
uses the actual admitted wrench bound. Only total work is constrained to be
dissipative: velocity weakening can transfer translational energy into spin.

Four focused native regressions cover analytic affine sliding/spin stopping and
travel, exact Coulomb reduction, reusable curve trials, spin gain with total
energy loss, unrelated-body preservation, geometry/normal refusal, cancellation
and deterministic retry. The four existing curved native and six exact native
trajectory tests also pass. Poses, load, radius, geometry and inertia remain held;
source/temperature history, normal evolution and external drive keep their
separate owners. Spatial pressure and physical-time errors require separate
refinement; noncoercive domains retain the constitutive refusal.

# CONTRACT: fs-mbd

## Curved native Hertz motion and reusable coefficient trials

`tangential_friction::hertz::transient` applies the original pressure-resolved
curved slip/spin trajectory to two finite native bodies. Actual masses, full
inertia tensors, contact arms and tangent frame supply the complete mobility.
The original physical-power quadrature and bounded extinction tail provide
heat; the existing compensated pair-impulse owner applies the integrated point
and free-couple impulses once. Actual native endpoint rates, body/pair momentum,
kinetic change and impulse work independently check the candidate before either
body publishes. An average wrench is not used as an endpoint friction law.

A conservative bound on the original normal-response row must preserve held
closed normal motion throughout the entire path, including the tail. Normal
evolution, external drive, changing geometry and pose integration remain outside
this frozen-contact branch. The numerical prefix estimates and explicit tail
time/travel bounds remain visible in the original constitutive report.

`PreparedRigidHertzTrajectory` and `PreparedRigidHertzTransient` capture the
original pair and full response once. Every coefficient trial starts from the
same accepted momenta and must retain the original normal load and radius.
Thermal/source identity stays with the enclosing source owner. Their contact
point getters expose the exact geometry used by native mechanics.

The native `hertz_transient` tests compare independent planar response and RK4
work/travel, analytical inverse-mobility stopping heat/impulse, stationary dwell,
tangent/world rotations, coefficient retries, unaffected bodies and atomic
budget/cancellation/normal-path refusal. The exact trajectory tests also check
prepared trials against fresh preparation and reject changed patch geometry.

## Finite-duration Hertz stopping trajectories

`tangential_friction::hertz::trajectory` stages or applies the original Hertz
shear/spin trajectory to two finite native bodies through an invariant-direction
stopping event and its stationary remainder. Actual masses, full inertia,
contact arms and tangent frame supply the complete mobility; one original
pressure-resolved kinetic-law evaluation supplies the moving phase. Only that
phase contributes impulse, displacement and physical loss. The existing
compensated pair-impulse owner applies each body's resultants once.

Actual native endpoint rates, moving-phase travel, point/couple work and
linear/angular momentum independently check the candidate before publication.
Normal rates remain exposed to the original normal owner. Poses, geometry, load
and material coefficients stay fixed; external drive and curved slip/spin
trajectories are outside this exact branch. Refusal or cancellation publishes
neither body's candidate. The `hertz_trajectory` tests cover analytic sliding
recoil, reduced-inertia spin heating, off-diagonal response, tangent-frame
covariance, split intervals, unaffected bodies and atomic refusal.

## Exact held isotropic drag

`linear_drag::IsotropicLinearDrag` advances the original native pose and momenta
under a fixed world force, uniform gravity, held translational bath velocity and
nonnegative isotropic translation/rotation resistance. The body must have exactly
equal principal moments and no applied couple. Translation and its displacement
follow the affine exponential solution; free spin decays about a fixed axis and
the attitude uses its integrated angular velocity. The shared `fs-math::decay`
owner retains tiny losses and representable late tails. There is no relaxation
number cutoff, midpoint sign reversal, speed clamp, inertia projection or new
material law; changing sources, lift, contact and general forced/tumbling
rotation remain outside this explicit reduction.

The receipt integrates drag impulses, held force/gravity work, maintaining-bath
work and translation/spin heat independently of represented-state energy.
Translation heat uses the nonnegative mean-square plus exponential-variance
integral, including its regular zero-drag limit. The shared native kinetic-change
owner and actual returned coordinates determine storage and numerical residuals;
lost coordinate or momentum increments remain visible rather than becoming heat.
`tests/linear_drag.rs` covers independent analytic driven motion/work, stiff
entrainment, semigroup composition, terminal settling, a tiny unresolved state
change, zero drag, Galilean covariance and unsupported-input refusal.

### Prescribed force ramps with simultaneous drag

`IsotropicLinearDrag::step_affine_force` uses the same initial-wrench/final-world-
force convention as the native integrator's affine-force entry. It integrates
the prescribed ramp and drag together, retaining gravity, held bath velocity,
fixed isotropic inertia and the original exact free-spin flow. A held force
delegates to the existing held entry. No separate force kick, terminal-speed
approximation or relaxation-number restriction replaces the actual path.

Applied work uses both the velocity integral and its first time moment; even
zero net applied impulse can transfer nonzero work. Independent slip heat uses
the positive squared projections of the affine/exponential trajectory onto
constant and linear time, plus its nonnegative remaining variance. The short-
interval moments have regular series with residual variance tending to 1/720;
the stiff branch uses the decaying-exponential basis. Zero translation drag
uses exact polynomial motion without reading bath velocity into the solution.
The original native pose, spin, kinetic-storage and impulse/work-residual
owners publish the result. Represented-state rounding remains a defect.

`tests/linear_drag/affine.rs` checks independent closed-form velocities and
power quadrature, a zero-impulse work transfer, unequal-interval composition,
small-exposure reversal heat, zero-drag ballistic work, Galilean source work,
held-entry parity and unsupported inputs. Endpoint force samples do not define
an arbitrary feedback law. A material consumer must admit interior speed
extrema; the held-force endpoint-only Reynolds argument does not extend to a
ramp. Accelerating-fluid pressure, added mass and history remain with their
existing physical owners; this entry changes the prescribed force only.

## Relative contact velocity

The shared articulated normal, tangential and bilateral velocity rows accumulate
their original signed wrench/twist scalar products with `fs_math::dd::Dd`,
reducing translation separately from rotation so the counterpart and prescribed
support contributions cancel before spin is added. Both
geometry preparation and the actual native endpoint query use this reduction.
Small rotational slip retained in the native twists therefore survives a much
larger common translation. Each row has fixed work and no allocation; original
frames, material witnesses and inverse-mass actions still determine the motion.

Focused G0/G1/G3 tests in `tests/articulated_contact_geometry.rs` check independent
rigid-point normal/tangent velocities, body and prescribed-plane counterparts,
ordinary boosts and world rotations, and an analytic off-center impact through
the native solver. An extreme exactly represented boost exercises cancellation:
if its required impulse endpoint is unrepresentable, the event refuses rather
than being misclassified as stationary. The reduction cannot reconstruct motion
already rounded out of an input twist or support projection, and does not alter
the existing endpoint admission tolerances. The public `relative_point_speed`
entry point validates consumed row/twist scalars and indices before delegating
to this same reducer; it does not certify geometry or pose correspondence.

The finite-body patch adapters use `PointKinematics::relative_velocity_world`
for initial slip, endpoint motion and conjugate work. This covers partial slip,
Hertz and adhesive shear/spin, rolling, directional friction and Mindlin contact,
including shared-body islands. Common COM translation cancels before either
rotational point contribution is added; subtracting rounded absolute point
velocities is not a relative contact query. The checked fixed-work reduction
propagates the original dynamics refusal before state publication and changes
neither the friction law nor the body momentum admission gates. The rotational
slip cases in `tests/tangential_friction.rs` and `tests/hertz_friction.rs` compare
rest-frame history, traction, recoil and heat under common tangential drift
using exactly representable impulses and an independent analytic Hertz loss.

`PointKinematics` also retains `relative_point_velocity_world`, the body-coordinate
rate rotated into world axes. Body-fixed material queries initialize it to zero;
`PreparedRigidPointMotion::relative_point_motion` records the prescribed rate
without changing the carrier COM or angular velocity. Both relative-velocity
queries include this component before their final rounding. The moving-point
absolute query uses the same reduction, so cancellation between translation and
prescribed motion does not erase a smaller rotational contribution. Fixed-point
queries retain their original arithmetic path. G1/G3 cases in
`tests/relative_point_motion.rs` check analytic rotated-coordinate derivatives,
relative motion between two points and a support, common-translation invariance,
and cancellation in the absolute query. A migrating geometric location remains
distinct from either body's material surface velocity in a contact law.

## Physical torque products

The native `Vec3::cross` compensates both floating-point product tails before
forming each moment. For finite product differences this preserves antisymmetry
and exact cancellation of parallel binary-scaled inputs, without a tolerance
clamp. Near-central impulses retain small resolved angular momentum instead of
losing it in two large products. The operation is allocation-free and fixed
work; nonfinite differences retain the previous fused-operation range behavior.
The G0/G1 cases in `tests/momentum_reductions.rs` cover native central impulses,
axial rigid-point motion, and all three nearly central torque components against
an independent exact integer oracle.

## Shared native free rotation

The explicit kick-drift, symplectic and staged Verlet rigid-body paths use
`fs-time::lie::gyrostat::InertiaTensor3::momentum_step` with zero internal
rotor momentum. The original admitted principal moments and fixed
principal-to-body frame enter that shared owner directly, without another
eigensolve or a momentum-to-velocity round trip. The common Hamiltonian
splitting advances attitude and body momentum together, including full
nonprincipal inertia and the exact axisymmetric-top special case.

The MBD boundary retains its normalized canonical quaternion and original
force, torque, translation and diagnostic conventions. The implicit
midpoint/Cayley path remains a separate integration method. The energy test
entry includes an independent continuous rotated symmetric-top reference for
all three explicit paths, with spatial momentum and kinetic-energy checks.
Native group admission and representability limits apply to the shared drift.

The native midpoint, Cayley, kick-drift and Verlet translation stages retain
the product tail of `mass * gravity` until it has combined with the applied
world force. That compensated resultant drives the momentum impulse and COM
drift, including separately evaluated endpoint forces. Nearly balanced loads
therefore retain their physical residual acceleration. The instantaneous
`prepare_point_motion` query retains the corresponding `force / mass` division
remainder until adding gravity; it does not require a representable weight
when the acceleration is representable. These fixed-size operations allocate
nothing and leave the independent force/gravity work audit intact. The
integer-oracle regression in `tests/momentum_reductions.rs` checks all native
translation lanes, signed axes, binary unit rescaling, varying endpoint loads
and derivative consistency. Representability limits of the returned `f64`
state remain; compensated arithmetic does not project momentum or energy.

## Fourth-order free rigid-body motion

`RigidBodyIntegrator::step_free_fourth_order` composes the shared compensated
principal-axis Hamiltonian flows with fourth-order symmetric coefficients.
Full nonprincipal inertia uses the cached original principal frame; the hot
path allocates nothing and performs no factorization or nonlinear solve.
Seven merged residual axis flows and one commuting spherical flow advance
asymmetric rotation. Spherical and axisymmetric bodies retain the original
exact flow. Every attitude increment transports angular momentum by its
inverse, preserving spatial spin to roundoff. The existing compensated native
translation remains exact for the integrator's configured uniform gravity.

This optional entry accepts no external wrench, contact or dissipative law;
its negative internal stages belong only to smooth free rotation. It does not
change existing step defaults or upgrade the order of a caller's surrounding
load splitting. Finite asymmetric energy and phase require timestep refinement;
no exact-energy or discrete-adjoint claim is made. The independent Euler-flow
order and equal-flow-count accuracy tests live at the shared `fs-time` owner.
Native `tests/angular_momentum_stages.rs` additionally checks full-inertia
linear/spin/orbital momentum, long-trajectory reversibility, analytic gravity
ballistics and the external orbital angular impulse from uniform gravity.


## Prescribed affine-force motion and work

`RigidBodyIntegrator::step_affine_force` accepts a world-force history specified
by its initial and final vectors, a held body torque, fixed physical inertia,
and the integrator's constant gravity. It analytically integrates the affine
force's momentum impulse and displacement first moment. For dimensionless time
`s=t/h`, the translated momentum is
`p(s)=p0+h*((s-s²/2)*F0+s²/2*F1+s*mass*g)`.
This is an explicit source history, not an exact solve for a state-dependent
force merely sampled at two endpoints.

Signed applied work independently integrates `F(s) . p(s)/mass` with the affine
force basis. Its gravitational cross term cannot be replaced by mean force
times net displacement: a force ramp can supply nonzero work even when both
its net impulse and the body's net displacement vanish. The gravity-work
account uses the retained path, while gravitational potential change uses the
represented endpoint positions. Kinetic storage likewise comes from the actual
endpoint momenta. Lost state increments and the original rotational splitting
error therefore remain visible in the mechanical work residual; none becomes
heat, a source correction or a state projection.

The original compensated held-body-torque kicks and full-inertia free drift
own rotation unchanged. Each step uses fixed stack storage and no allocation,
factorization or additional nonlinear solve. Callers may compose adjacent
affine segments at their shared source endpoint. Arbitrary force histories,
force discontinuity discovery, changing inertia, contacts and exact asymmetric
rotation are not claimed. `tests/driven_work.rs` covers analytic 3-D forced
motion, an independent Simpson power integral, zero-impulse energy transfer,
braking, split continuation, full-inertia rotational parity and its unchanged
defect, retained source work beside rounded storage, and input refusals.
This advances `frankensim-ext-mbd-core-vqqt`; that broader bead remains open.

## Simultaneous internal-wrench networks and finite motion

`spatial_loads::internal_impulse::network::apply_internal_wrench_network`
extends the original internal connector to cycles and shared finite bodies.
Each connector uses the original attachment wrench transport, then all linear
and central angular impulses accumulate in a compensated bank before each
original native body updates once. The original per-body representability
audit and whole-network linear, spin-plus-orbital and work gates remain active.
Opposite undelivered body resultants cannot hide behind zero net momentum.

Every connector's signed work uses the complete simultaneous before/after
midpoint velocities. This retains kinetic cross terms when multiple sources
act on one body, with no arbitrary sequential source attribution. A second
linear connector pass reuses original geometry; the body bank becomes midpoint
kinematics. Work and owned payload are O(bodies + connectors), with body,
connector and complete peak-byte limits admitted before fallible allocation.

`spatial_loads::internal_motion::network::step_internal_wrench_network` wraps
one original zero-load free step per body between two complete network stages.
The transported couples follow actual endpoint attachment geometry, and all
connectors and isolated bodies share one clock. Each free body's actual drift
is admitted separately before whole-step momentum and source-work admission.
Peak storage includes both receipts and the active native workspace. Any
refusal or observed cancellation publishes no partial bodies or work.

These are prescribed world-held internal sources with explicit supplied or
recovered mechanical work. They do not infer joint/contact laws, elastic
storage, motor efficiency or heat. The symmetric finite trajectory and source
work require step refinement for general rotating attachments.

Five G1/G3/G4 groups in `tests/internal_wrench_network.rs` independently check
shared-body kinetic cross terms, full-tensor cyclic wrench algebra, body
relabeling and frame/common-boost covariance, exact unequal-mass force/couple
chains, and quadratic convergence to a simultaneous three-body Newton/Euler
ODE with individual source-work quadratures. They also cover reversal,
single-pair and zero-source native parity, complete payload/count bounds,
swallowed resultants, late refusal and every observed cancellation point.
All five groups and eleven original pair impulse/motion groups pass in the
bounded actual-source native lane; full Cargo/DSR is not claimed.

## Finite prescribed internal-wrench motion

`spatial_loads::internal_motion::step_internal_wrench` advances two original
finite bodies through a held world force/free couple using native internal
half-impulses around the original zero-load rigid-body drift. Both actual
endpoint attachment geometries determine the transported opposite couple.
Original full inertia, quaternion motion and implicit Euler-top solve remain
authoritative; the bounded operation uses two impulse evaluations and two
free-body steps without allocation or hidden substeps.

This symmetric second-order trajectory separately gates each free body's
actual momentum and energy drift, then the pair's whole-step linear and
spin-plus-orbital change and actual kinetic change against signed actuator work.
Double-double endpoint reductions prevent opposite individual drift defects
from canceling. The returned work is a discrete conjugate half-impulse
quadrature; free-drift energy defects remain visible in the work gate, never
heat. The actuator may supply or recover energy. No spring/contact law,
passivity or symplecticity is implied. Invalid duration, unrepresentable
nonzero half-impulses, native solve refusal, budget failure or cancellation
returns no partial endpoint/account.

Six G1/G3/G4 groups pass in the real source-selected native lane: unequal-mass
constant-force motion, finite motor-angle exchange/regeneration, full-inertia
off-centre convergence against an independent refined continuous RK4 reference,
conservation about multiple origins, rotated/translated/common-boost covariance
and reversal, exact zero-drive reuse with independent body drift gates, and
every cancellation/representability refusal. Orientation, spin and work errors
reduce approximately fourfold under each step halving. Rustfmt check also
passes; full Cargo/DSR is not claimed.

## Loaded finite internal actuation

`spatial_loads::internal_motion::step_internal_wrench_with_loads` and
`internal_motion::network::step_internal_wrench_network_with_loads` compose the
simultaneous internal half-kicks with one loaded original rigid-body midpoint
step per body. `HeldExternalLoads` borrows one world-held COM force/body-held
free couple per original body plus uniform gravity. Isolated bodies still
receive the environment; mismatched load banks refuse. Original masses, full
inertias, poses and nonlinear rotation remain authoritative.

With no body couple the original default native step retains its exact
exponential special cases. A nonzero body-held external couple selects the
existing Cayley owner, giving independent spin impulse
`dt*(R0+R1)/2*tau_body`; external orbital impulse is
`(x_mid-origin) cross dt*(F_world+m*g)`. Products, source reductions and
original-state momentum differences use compensated arithmetic. Source
impulses are never inferred from observed momentum changes.

Native `StepWorkDiagnostics` supplies force work on actual COM displacement,
conjugate body-couple work and gravitational work. These remain separate from
each internal connector's signed actuator work, including regeneration. Each
body's middle-step source balance must close before the complete pair/network
gate admits `delta(K+V_g)-W_external-W_actuator`. Numerical drift defects remain
signed and are never deposited as heat. Internal kicks do not move poses, so
`delta V_g=-W_g` covers the complete trajectory.

The composition is second order for smooth held commands; rotations and
continuous actuator work require timestep refinement, especially near the
Cayley half-turn. Pair motion remains allocation-free, and network work/storage
remains linear in bodies plus links with the existing complete-peak budgets and
atomic refusal. Unloaded entry points retain their native motion. This does not
infer constraints, external point-force arms, feedback control or a thermal law.

Focused tests cover closed-form unequal-mass recoil/free fall, loaded rotor
generation, shared-body force chains and isolated bodies, full-tensor
trajectory/source-work refinement against the independent Newton/Euler ODE,
invalid sources, load counts and zero-load parity. These added tests have not
been executed in the current resource-limited environment; the broad rigid-body
lane remains open.

### Spatial point loads on finite actuator networks

Finite internal actuation accepts `HeldSpatialLoads` through
`step_internal_wrench_with_spatial_loads` and
`network::step_internal_wrench_network_with_spatial_loads`. One borrowed
`PreparedSpatialLoads` per original body retains world-held COM/material-point
forces plus body/world free couples. Existing simultaneous actuator half-kicks
surround the original bounded full-inertia spatial midpoint owner. Native point
geometry, nonlinear replay and Cayley limits remain authoritative. Exact zero
prepared loads retain the original free-body chart.

Source force work is the actual material-point force-displacement contraction,
including rotation of the lever arm. Free-couple, gravitational and internal
actuator work remain distinct. Source angular impulse uses original compensated
prepared first moments and endpoint rotation averages about the declared
origin, independently of observed momentum change. Each native middle motion
and the complete network must pass the original dimensional work/momentum
gates. Embedded native step receipts use the finite-step torque pullback;
`external_loads` carries the physical force/free-couple work classification.

Pair motion remains allocation-free. Network work/storage remains linear in
bodies plus connectors, independent of original point count after preparation.
Continuous phase/source work require time refinement. No follower force,
contact, passive connector or callback law is inferred. Focused tests add
independent endpoint point work and quaternion-derived couple work, source
impulses and origin covariance, full-tensor three-body trajectory/power
refinement, native limit refusal, cancellation and zero-load parity. These
added native tests have not been executed in the current resource-limited
environment.

## Native internal spatial-wrench impulse transfer

`spatial_loads::internal_impulse::apply_internal_wrench_impulse` transfers a
declared force/couple impulse between two original finite bodies at their actual
material attachments. The second attachment receives `(J, C_B)`; the first
receives `(-J, -C_B - (x_B-x_A) cross J)`. The arm term is the physical moment
transport required for an internal connector between noncoincident points. It
is applied before the original full-inertia momentum update; no state is
projected and the existing independent point-impulse semantics remain intact.

Each body retains its own native mass, full tensor, pose and momentum frames.
Explicit dimensional gates check each commanded versus representable impulse,
the pair's actual total linear and orbital-plus-spin angular change about a
declared origin, and actual quadratic kinetic change against conjugate connector
work. Constant-size double-double reductions retain the original transported
central wrench, native frame conversion and actual endpoint differences. They
avoid subtracting large rounded background momenta or energies. Opposite
swallowed kicks, erased old momenta, and lost small free couples therefore
refuse even when their net balance would appear to close. No heap storage is used.

The connector can supply or recover mechanical energy. Its signed impulse work
does not imply passivity, a motor constitutive law, frictional heat, a contact
response or a finite-duration joint solve. Cancellation, nonfinite arithmetic
and unresolved physical budgets return no changed body or partial account.

Five focused G1/G3/G4 groups cover noncoincident full-inertia recoil and angular
conservation about different origins, the coincident point-impulse limit, pure
motor-couple spin exchange and energy recovery, rotated/translated/common-boost
frames, a resolved 1 J transfer against large background energy, refusal of two
swallowed 1 N s kicks, erased old linear/spin momentum, lost small free couples
behind large attachment moments, frame-conversion losses, and every cancellation
phase with deterministic retry.
They pass in the bounded native source-selected lane against the unchanged
original momentum owner. This is not a Cargo/DSR workspace validation claim.

## Supported finite rotational and translational attachments

The supported spatial-motion entries also accept the original `SpatialBushing`
law. Its anisotropic translational and principal-log rotational strains,
SO(3) force/couple transformation, stiffness, damping, coefficient admission
and chart boundary remain owned by the existing bushing module. The same
initial/midpoint geometry and effort evaluations feed supported dynamics;
both axial connections and bushings contribute to the returned spring storage,
internal work and damping accounts.

The endpoint must remain inside the native active rotational chart. Rejected
chart crossings or cancellation publish no partial motion. Empty bushing
slices retain the ordinary supported spatial path. This extension uses the
explicit supported stepper and does not introduce stiff or implicit bushing
integration. Its focused regression entry is
`cargo test -p fs-mbd --test articulated_prescribed_bushings`.

## Shared spatial loads on prescribed supports

`articulated::motion::spatial_loads::prescribed` connects the original
`SpatialMotionLoads`, world point forces/free couples, and finite axial
Kelvin–Voigt connections to the supported motion owner. Both load and effort
evaluation use the actual initial and midpoint states. The shared geometry
and constitutive evaluator are reused; no alternative force law or large
surrogate support mass is introduced.

All connection wrenches enter native ABA and RNEA, including their effect on
support reaction. External work, actuator work, signed support work, internal
connection work, actual spring-storage change and nonnegative midpoint damping
loss have separate accounts. The combined residual subtracts support work
exactly once. Damping is not deposited into a thermal model implicitly.
Final connection geometry and the last cancellation poll must also succeed
before the endpoint is returned.

The existing explicit timestep and connection-mode resolution requirements
apply. Rotating supports preserve the world direction and actual application
point of external loads; scalar root joints remain admissible. The targeted
regression entry is `cargo test -p fs-mbd --test articulated_prescribed_spatial`.

## Local response about supported mechanical equilibria

`articulated::motion::spatial_loads::prescribed::equilibrium::linearized`
prepares local mechanical response about an original stationary or uniformly
rotating supported equilibrium. Preparation independently reevaluates the
physical force balance against the immutable model, loads and material laws;
a public root report cannot replace that admission. Constant held actuator
efforts, native axial connections and diagonal/coupled bushings retain their
original meaning. The rotating owner's world-load symmetry checks are shared.

The original joint coordinates obey
`M*delta_qdd + C*delta_qd + K*delta_q = delta_effort`.
Matrix-free inertia uses native inverse dynamics with the direction in the
acceleration argument and all independent bias terms zero. Position and rate
actions differentiate the same original required-minus-applied effort balance.
Joint accelerations stay zero while joint velocities are perturbed, preserving
the original velocity bias without confusing it with inertia. K includes
geometric preload, follower loading and centrifugal softening. C includes
material damping and Coriolis/gyroscopic terms. Neither operator is projected
onto symmetry, positive stiffness or proportional damping, and an unstable
root remains unstable.

The shared bounded central directional-action owner compares successive
stencils in explicit position/rate and effort scales. Physical coordinates are
anchored at the admitted state and incremented with fused multiply-add.
Compensated requested-versus-represented offset checks enforce an explicit
relative stencil-distortion allowance. Overflow, nonzero values erased by
scaling, original physical-domain refusals and unresolved stencils return
without a partial action. Reported stencil disagreement is an estimate, not a
certified derivative error, smoothness proof or branch-stability certificate.
The original state, clock, loads and material history are never advanced.

An explicitly bounded dense export supports at most 64 joint coordinates.
Its mass block comes directly from the existing generalized composite-rigid-body
owner, preserving that owner's symmetric assembly rather than repairing
independently sampled mass columns. The export retains full K and C, physical
stencil diagnostics and an explicit source-evaluation budget. Existing fs-modal
quadratic pole and prepared-mode owners consume the same M/C/K for poles,
mode shapes, residues and sensitivities; no second eigensolver is introduced.
Larger admitted trees can use matrix-free actions within original tree limits.

Analytical integration cases cover a damped spinning radial slider, transverse
preload stiffness, inverted-pendulum growth, a two-coordinate Coriolis coupling,
original balance/budget/cancellation refusals and coordinate quantization.
These new integration cases remain unexecuted because the execution workspace
disconnected. The shared derivative kernel's six native cases passed before
that outage; its final underflow hardening and this full mechanical integration
were reviewed afterward. Arbitrary effort callbacks, contact switches,
closed-loop constraints and nonuniformly rotating operating points are outside
this local supported-tree response API.

## Relative equilibrium under supported uniform rotation

`articulated::motion::spatial_loads::prescribed::equilibrium::rotating::equilibrate_uniform_rotation`
solves a steady joint shape in the original uniformly rotating supported tree.
It reuses the existing scaled Newton--Krylov force-balance owner and the exact
same native zero-joint-acceleration inverse-dynamics residual. The actual base
angular velocity remains present in kinematics and RNEA, so centrifugal and
gyroscopic effects arise from the original mass properties. No separate
centrifugal force, rotating potential or dense mass/tangent matrix is added.

The base has a fixed support origin, zero linear velocity and zero spatial
acceleration. Joint rates are zero. Internal axial connections and diagonal or
coupled bushings retain their original moving attachment geometry and material
laws. Held body wrenches and generalized efforts co-rotate. World gravity,
point forces and couples must be invariant about the original rotation axis;
scaled compensated cross products admit only machine-roundoff axis error and
never project or modify the applied vectors. Arbitrary actuator callbacks are
outside this symmetry-admitted entry. Zero spin uses the stationary mode.

The returned state retains the original support motion and clock and the
independently checked physical joint residuals. The support reaction is still
the original root wrench. Its compensated dot product with the actual base
twist reports drive power in watts into the tree; this may be nonzero under a
co-rotating drag. Attachment storage remains the original constitutive energy.
Neither drive power nor displacement from the seed invents a heat destination
or a relaxation-work history.

`UniformRotationFrame` supplies `R(t)=R(0) Exp(omega_body*t)` at the original
represented origin. `base_at_elapsed` supplies native midpoint and endpoint
boundary samples to the existing prescribed-motion evolution owner, without
accumulating orientation steps. Original SE(3) representation and geometric
admission remain authoritative: translation can round and a distant rotated
pose can be refused. Overflowed or completely underflowed nonzero rotation
increments refuse. Zero elapsed time and zero spin retain the seed pose.

Physical row scales, joint limits, finite source/material charts, bounded
residual/Krylov work and cancellation remain the same as stationary equilibrium.
Refusal leaves the seed unchanged. This is a driven relative equilibrium, not
a zero-kinetic-energy state, a free-body angular-momentum solve, a stability
certificate, a contact/sticking solver or a closed-loop mechanism solution.

Four native G0/G1 production-frame checks passed against the actual fs-ga
implementation: arbitrary-axis trajectory derivatives, invariant load scaling,
zero-motion/invalid-time behavior and retention of original distant-pose
refusals. Five authored integration regressions cover an independently solved
centrifugal spring extension, original prescribed-motion continuation and drive
power, a tilted axis, an asymmetric rotor's gyroscopic reaction, stationary
parity, symmetry refusals and an unsupported supercritical spring branch.
The full articulated integration target has not run in this environment because
the genuine asupersync dependency exceeds the available compilation memory.
The small native frame check is not a whole-MBD or whole-workspace validation.

## Static equilibrium on supported articulated mechanisms

`articulated::motion::spatial_loads::prescribed::equilibrium` solves the
original zero-rate RNEA holding effort minus actual actuator effort. World
gravity, link-frame wrenches, world point forces/couples, finite axial links
and diagonal or full passive six-axis bushings use the same geometry and
constitutive evaluator as supported time evolution. `equilibrate_with_effort`
accepts the original smooth supported-motion effort callback; its elapsed time
is zero and its absolute clock is held. Repeated trial calls must not advance
accepted material or control history.

Every joint has explicit positive coordinate and effort scales in its original
SI units. Newton--Krylov and directional finite differences belong to
`fs-solver`; there is no assembled generalized mass or Jacobian. Initial
coordinates are preserved exactly by solving scaled increments. Independent
final RNEA evaluation must meet the absolute scaled force 2-norm tolerance;
relative initial-residual stopping cannot weaken it. The result retains each
physical joint residual, required/applied efforts, original attachment storage
and the complete support wrench. Support extraction is shared with the original
time-evolution owner, including scalar-root actuator transmission.

The seed must have zero rates and base acceleration. Returned coordinates and
zero rates explicitly initialize the native state while preserving its clock
and support pose. No elapsed time, relaxation work, damping heat or loss is
inferred from movement between the seed and solution. The result can continue
through the unchanged supported time integrators. Native joint limits and
attachment chart refusals apply to every trial and numerical stencil.

Work is bounded by 4096 links, the original attachment inventory bounds,
64 Newton attempts, 64 Krylov vectors and 16384 residual evaluations. Every
evaluation, including final verification, consumes the caller's budget and
polls cancellation around bounded native passes. All-fixed models recover
support with no artificial zero-dimensional solve. A loaded singular mechanism
does not publish an equilibrium; a balanced neutral or unstable root may be
accepted without a stability, uniqueness or reachability claim. Contact,
hard-stop reactions, Coulomb sticking, moving-support relative equilibria and
closed-loop constraints remain separate owners. Numerical derivatives are not
analytic tangents or adjoints.

`cargo test -p fs-mbd --test articulated_equilibrium` covers G1 loaded axial
extension and support force, finite-angle gravitational torque, reciprocal
translation/rotation coupling and physical rescaling, plus G0 no-drift native
continuation, original effort callbacks, all-fixed reactions, singular-load
refusal, budget exhaustion, cancellation and seed immutability.

## Continuous prescribed-base articulated motion

`articulated::motion::prescribed` advances the original compact joint state
with a fixed or externally driven support. `PrescribedBaseMotionState` retains
physical time and the complete native `BaseState`; `PrescribedBasePath` supplies
the midpoint and endpoint of the same smooth support trajectory. The original
Featherstone spatial-acceleration convention applies. Pose/twist/acceleration
derivative consistency is an input obligation, not a certified interpolation.

Both `step_prescribed_base` and `step_prescribed_base_with_effort` use native
linear-time ABA for the explicit midpoint stages. Native FK and midpoint RNEA
recover the complete support force and couple, including transmission through a
scalar root actuator. The root transform comes from the original local joint
geometry. Support work pairs this wrench with base twist; actuator work pairs
effort with relative joint rate. The shared held-load carrier retains its
original link-frame wrench and world-gravity conventions.

Endpoint storage, polarized kinetic-energy change, exact uniform-gravity
displacement work, signed support work, and actuator work remain separate.
The kinetic-work defect is numerical/path error, not physical dissipation.
Joint position/rate/effort admissions and bounded cancellation apply before
publication. Refusal preserves borrowed native state; callback-owned side
effects cannot be rolled back. The accepted rates from
`apply_prescribed_base_impulses` continue directly without converting an impulse
into a finite force or resetting configuration.

This is a smooth explicit second-order lane, requiring time resolution of
its fastest dynamics. It is not a stiff, symplectic, energy-exact, unilateral
contact or closed-loop constraint integrator. The targeted regression entry is
`cargo test -p fs-mbd --test articulated_prescribed_motion`; it covers analytic
support/actuator work, driven supports, native impulse continuation and refusal.
This increment advances `frankensim-ext-mbd-core-vqqt`; it does not close the
bead's broader integrator program.

## Momentum-consistent fixed-plane candidate projection

`rigid_rattle::PointPlaneRattle::project_candidate` carries world linear
momentum and world central spin through its configuration repair before
applying the original full-inertia normal point impulse. COM position is
corrected strictly along the fixed world normal, while body attitude uses
its own exponential. A changed body basis must not rotate the stored world
momentum, and a finite rotational repair must not displace the COM tangentially.
The composed projection preserves plane-tangent linear momentum and total
angular momentum along the plane normal to roundoff. Normal COM relocation
still changes the other orbital components; it does not manufacture spin to
cancel them. Anisotropic kinetic energy can change during configuration repair.
This is a bilateral candidate projector, not a full variational RATTLE time
integrator, contact law, or an energy-conservation claim. Configuration, tangent,
and composed projection failures leave their caller-owned states unchanged.
Focused G0/G3/G4 tests cover full-inertia point response, plane symmetries,
world rotations, COM-only repair, and nonlinear/endpoint refusal.

## Physical acceleration at articulated material points

`InverseDynamics::prepare_point_motion` and `ForwardDynamics::prepare_point_motion`
use the original matching `Kinematics` and original gravity to restore physical
acceleration from prescribed-base RNEA/ABA's gravity-shifted spatial receipt.
`FreeFloatingForwardDynamics::prepare_point_motion` accepts its already physical
spatial acceleration without adding gravity again. Both apply the necessary
`omega cross v` conversion from spatial to Cartesian origin acceleration.
The prepared link query retains `alpha cross r` and `omega cross (omega cross r)`
for each original body-fixed material point. Native joint Coriolis/centrifugal
effects remain in the original solver's acceleration, without another tree or
inertia solve. Preparation and queries are fixed-work and allocation-free.

`ArticulatedPointMotion::specific_acceleration_world` separately returns ideal
accelerometer specific force `a - g`: a stationary supported point reads `-g`,
whereas a nonrotating point in free fall reads zero. No impact acceleration,
sensor dynamics, sliding attachment, missing reaction solve or time integration
is inferred. Callers must supply receipts from the same model, state and instant;
shape/index and finite-arithmetic admission cannot prove that provenance.
G1/G3 tests cover an analytic rotating radial slider, independent differentiation
of native forward kinematics, rotated frames, common velocity boosts, gravity
conventions and malformed/unrepresentable motion refusal.

## Prescribed-base articulated impulse response

`articulated::apply_prescribed_base_impulses` applies body and joint impulses
with the original base pose and twist held fixed. The existing ABA recursion
solves the constrained joint inverse mass in O(n) work/storage; its condensed
six-dimensional residual supplies the support impulse on the retained tree.
No heavy free-base approximation, timestep scaling or velocity projection is
used. Scalar root joints and fully fixed trees are admitted. Position and actual
endpoint speed limits apply; instantaneous impulses are never compared with
force/torque effort limits. Smooth acceleration and gravity contribute no
finite impulse over this zero-duration event. A nonzero solved joint jump that
rounds away completely refuses, including a changed joint rate erased by an
unchanged large parent twist; an ideal constrained zero jump remains valid.

The receipt retains the reaction about the base origin in base coordinates and
about the world origin in world coordinates, including off-axis couples and
root-actuator transmission. Actual native link velocity increments independently
give linear/angular momentum changes and polarized kinetic change. External
body and relative joint work remain separate from signed prescribed-support
work. A stationary support does zero work; a driven support can supply or absorb
energy. Residuals are reported without projecting or repairing the endpoint.

`articulated_impact::resolve_prescribed_base_plane_impact` uses the same scalar
Newton law and point kinematics as free-base impact, with the supported inverse
mass and support work. It admits one supplied contact configuration with
positive normal mobility; locked directions refuse. Both base boundaries check
the actual returned normal restitution on the relative impact-speed scale;
large common motion cannot admit a collision lost in endpoint rounding.
It does not detect impact time, advance pose, add compliance or allocate
dissipation to heat. The separate simultaneous owner below combines supported
and free trees. Native G0/G1/G3 tests use full-tensor off-centre
hinges, a coupled slider mass oracle, driven supports, scalar root actuation,
frame changes, a fully fixed tree, impulse units and actual endpoint limits.

The simultaneous free-tree normal adapter also consumes the contact owner's
unique-velocity result for redundant rows. Its individual impulses retain only
representative status. The contact owner's active preimpact velocity range
admission and native mechanical Gram structure determine total support work;
unresolved work refuses. Actual endpoint, momentum and energy audits remain
authoritative.

### Mixed free and prescribed articulated contact islands

`articulated_contacts::mixed::resolve_mixed_articulated_contacts`,
`resolve_mixed_articulated_frictional_contacts` and
`resolve_mixed_articulated_elliptic_contacts`
resolve simultaneous normal, circular-Coulomb and directional-Coulomb events
over one common set of free and prescribed-base trees. Each input declares
`ArticulatedContactBoundary`; the returned `MixedArticulatedImpulseResponse`
retains the corresponding native response in input-system order. Original
free-only entry points and their concrete response types remain available.

Each contact row probes only its involved trees using their original free or
prescribed ABA inverse mass. The full contact-space mobility retains shared-link,
joint/base and normal/tangent coupling. One law solve precedes one accumulated
physical impulse application per tree. A supported tree's base pose and twist
stay fixed while its unconstrained joints respond. A fully fixed tree can
participate in an impact against a mobile counterpart; an approaching contact
whose complete relative direction is locked refuses. No surrogate heavy base,
independent link mass or sequential contact clipping is substituted.

The global receipts include the full six-component prescribed-base reactions
in linear and world-origin angular momentum balance. `plane_work` and
`support_work` retain signed source contributions;
`boundary_work` is their sum. Driven translation and rotation can supply or
absorb energy. Contact dissipation remains distinct from this external work.
Under a unique-velocity law, the native Gram mobility and the contact owner's
work-range witness determine motion and total boundary work; individual contact
impulses, base reactions and plane/base work shares remain representatives.
Compatible rotating supports and independent contact islands are admitted.
The circular and directional owners also admit distinct harmless tangent belt
velocities when the corresponding friction coefficients are exactly zero;
ambiguous total relative work refuses in the original contact owner.

Actual returned link kinematics must reproduce the contact law on the relative
speed scale, independently of a large common support velocity. Original joint
endpoint limits, momentum/work audits and cancellation all precede publication
of the complete candidate. Configuration and borrowed inputs are unchanged.
Seven G0/G1/G3 groups cover full-tensor hinge reactions, driven support work,
redundant and permuted contacts, a supported-wall/free-body impact chain,
finite-mass projectile/hinge exchange, circular and directional sticking/sliding,
rotating support work, late limits/cancellation and unresolved stored speeds.
These are supplied-configuration zero-duration events, without collision
discovery, support compliance, time advancement or automatic heat allocation.

### Signed spatial momentum and impulse work reductions

Native rigid and articulated momentum receipts accumulate signed linear and
world-origin angular contributions with double-double reductions before the
final `f64` result. Spin, COM offset and world-origin orbital terms remain
separate through accumulation, preserving finite resultants and couples amid
large opposing momenta. Actual before/after velocity increments are paired for
linear/angular changes, and impulse work retains separate endpoint pairings
instead of first rounding a midpoint velocity. The same signed reduction
preserves separated-contact couples. This improves the numerical receipt of
the original stored `f64` state; it does not create extra physical precision,
project momentum or change the native inverse mass or impact law.

## Instantaneous rigid-point motion

`RigidBodyIntegrator::prepare_point_motion` captures one native smooth state,
constant mass/full inertia and instantaneous external wrench, adding configured
uniform gravity exactly once. The wrench excludes that gravity and must include
all other actual loads/reactions. Euler's original full-inertia action retains
the `omega cross L` term; `PreparedRigidPointMotion::point_motion` supplies native
position/velocity and tangential plus centripetal acceleration at each body-fixed
arm. Preparation shares the angular work across attachments; queries use fixed
work and no allocation, mutate no state and refuse nonfinite derived vectors.
They do not integrate, solve missing reactions, differentiate impacts, admit
changing inertia or describe sliding/deforming attachments. Acoustic coupling
may observe this motion but must supply its separate physical force on fluid;
gravity or the total mechanical wrench does not create an acoustic source.

`RigidPointMotion::specific_acceleration_world` uses the same ideal
accelerometer convention as articulated points: physical `a_point - gravity`,
in world m/s². Supported rest reads minus gravity and nonrotating free fall
reads zero. Gravity remains explicit even when a caller's external wrench
already includes weight. Sensor-axis rotation and measurement dynamics belong
to the consumer; `fs-flyer` uses this native observation in its six-DOF cue path.

## Exact linearly weighted tetrahedral inertia

`SpatialInertia::from_weighted_tetrahedron` accepts a positive total mass,
four body-frame vertices and four positive relative density weights. Barycentric
linear interpolation defines the mass distribution; common weight scaling
does not change it. Analytic degree-three moments supply the density-weighted
centre and full centroidal tensor. Nine centred positive dyads avoid cancellation
between large origin moments, with fixed work and no quadrature allocation.
The supplied mass is retained exactly, allowing affine mesh deformation to
transport reference matter without re-querying density or drifting mass.

Nonfinite/collapsed tetrahedra refuse even when the supplied mass is positive.
The original inertia owner admits the complete result. This numerical helper
has no material/source dependency, no mesh-overlap authority and no uncertainty
certificate. Source-density resolution and per-cell Jacobian transport live in
`fs-mbd-material::density_field`; its focused tests check independent Cartesian
integrals, all tetrahedron vertex orders, refinement, deformation and actual
native impulse/buoyancy response.

## Simultaneous resultant momentum impulses

`RigidBodyState::apply_momentum_impulse` accepts one world linear impulse and
one body-frame central angular impulse at fixed native pose. Original mass and
full inertia remain authoritative. Assemblers can sum all incident contact
impulses before this single update, avoiding contact-order-dependent sequential
rounding. `apply_impulse_at_body_point` delegates its native momentum update to
this primitive while retaining its actual point-work receipt.

The resultant receipt distinguishes commanded impulses from representable
momentum increments and computes kinetic change by polarization of those actual
increments. It reports linear, angular and work residuals without accepting or
repairing them; the composing contact owner supplies physical gates. This is
neither a contact law nor an angular-velocity state replacement.

The original simultaneous rigid-contact owner now uses this primitive for both
the circular and elliptic Coulomb paths. It accumulates every incident point
impulse and its ORIGINAL body-frame moment in canonical contact order before
one state update per body. Compensated reductions use the existing `fs-math`
double-double owner; no point is relocated and no force couple is discarded.
Native actual increments supply the total linear balance and
`x cross delta_p + R delta_L` angular balance about the declared world origin.
`kinetic_change_j` comes from polarized native increments, rather than
subtracting two large kinetic-energy totals. The separate total energies remain
diagnostic outputs. Unchanged large common translation or orbital momentum can
therefore no longer mask a small contact loss or cause a false balance refusal.
All original physical feasibility, endpoint-velocity, momentum and energy gates
retain their tolerances; this does not forgive an unrepresentable physical
impulse or infer a new contact law. G0/G3 regressions cover opposing simultaneous
impulses, canonical permutation, a common transverse boost, origin shifts and
normal restitution with unchanged controls.

The nonsmooth finite-patch tangential owner and the shared Hertz/adhesive
shear-spin island now use the same native resultant primitive. They sum every
point impulse and central couple before one body update, with compensated
component reductions. Their public `kinetic_change_j` and momentum balances
come from actual native increments; `kinetic_before_j` and `kinetic_after_j`
remain absolute diagnostics. Heat couplers must consume the change field,
since subtracting the diagnostics can discard friction losses under a common
translation. Physical heat, reversible shear storage, coefficient work and
numerical dissipation retain their separate original accounts and tolerances.

Focused G0/G3 checks retain a small original momentum and a small resultant
between opposing large contact impulses, and compare the original tangential
and shared Hertz solves with a common transverse/normal velocity. These are
fixed-geometry impulse phases, not new normal-contact, pose-integration or
constitutive laws; refusal still publishes no body or friction-history update.

Individual Hertz and source-backed adhesive pairs also expose the native
`kinetic_change_j` and use actual increments for their work and momentum
balances. Their fixed two-body resultant update stays on the stack, combining
the original point moment and central couple before rounding native momentum.
Adhesive temperature refresh closes represented storage with this kinetic
increment plus the normal-potential change; normal parameter work remains
separate from junction heat. Common-normal-velocity regressions cover Hertz
sliding, sourced adhesive spin and the real-source normal thermal refresh with
unchanged physical controls.

The shared articulated contact event owner also accumulates all incident
point-force and moment covectors with the existing `fs-math` double-double
arithmetic before rounding one spatial impulse per touched link. Sparse scratch
is bounded by the event's touched links; the original articulated inverse-mass
operator, contact geometry, constitutive laws and native endpoint update remain
authoritative. Signed prescribed-support work and actual native momentum/kinetic
increments retain compensation through balance subtraction. Opposing large
transfers therefore retain small physical resultants and support work
independently of contact/system order or unrelated common translation. Existing
feasibility and balance tolerances remain unchanged, and cancellation before
final admission publishes no endpoint. Analytic regressions specify force/torque
cancellation, boosted native balances, moving-support work cancellation through
a real normal-contact solve, and early/late transactional refusal.

## Linear rotor-bearing operating points

The public `fs_mbd::rotordynamics` module exposes the original operating-point,
rigid-rotor, transient and periodic owners described here. Finite internal shaft actuation
and relative rotor phase are available separately in
`fs_time::lie::gyrostat::rotors`; that fixed-axis angular model does not change
this bearing/bending lane's constant-speed physical assumptions.

`rotordynamics::RotorBearingModel` admits 1–32 original generalized coordinates
with full mass M, bearing velocity operator C, tangent stiffness K and the
skew gyroscopic coefficient G. M must be symmetric and numerically SPD; C/K
may be nonsymmetric or active and remain unchanged. A signed frozen spin rate
assembles `C + Omega G`. Coordinates and their conjugate forces/torques must
share the caller's declared reference frame and units, so `f dot q_dot` is power.

`RotorOperatingPoint::poles` delegates the complete quadratic pencil to the
existing scaled `fs-modal` owner; cancellation reaches its native complex QR
through `quadratic_eigenvalues_with_cancel`, with 120 sweeps per active window.
`campbell_samples` returns complete pole sets at 1–256 caller-supplied speeds
for those fixed coefficients. All complex poles, including active/unstable
ones, remain visible. `spectral_abscissa` retains the maximum computed pole
real part without an engineering deadband. `stability_transition_brackets`
uses only sign changes of that complete spectral abscissa and bounded bisection;
it never follows a sorted pole index, so active cross-coupled bearings can expose
a numerical onset interval even when QR order changes. Tangential touches with
no sampled sign change are intentionally not claimed.

`synchronous_critical_speed_brackets` likewise avoids manufactured mode
identity. At each nonnegative shaft speed it counts positive-frequency poles
below the declared synchronous line `omega_exc = order * Omega`. A change in
that integer count brackets one or more Campbell crossings. The routine bisects
the first count change and continues through the remainder of the supplied
interval, retaining unresolved simultaneous/nearby crossings as
`crossing_count > 1`. A changed positive-frequency pole inventory during
refinement refuses rather than reclassifying an overdamped/real-pole event as a
critical speed. This count invariant is insensitive to ordinary branch swaps,
but an even number of crossings can cancel between coarse sample endpoints;
the caller must resolve the speed grid.

These are deterministic **numerical brackets**, not certified eigenvalue or
coefficient-uncertainty enclosures. There are still no eigenvectors or persistent
physical mode identities, and no forced-response peak, flexible-shaft continuum
or Floquet certificate is inferred. Speed-dependent bearing coefficients require
a newly supplied physical operating point, not just another sample speed.

`rotordynamics::resonance` supplies the separate forced-response path.
`harmonic_response_with_derivative` differentiates the original full
`K - omega^2 M + i omega (C + Omega G)` equilibrium, retaining independent
forcing-frequency, shaft-speed and complex-force derivatives. The native
row-equilibrated LU is reused for the tangent; the original physical force
rows independently check both solves. Velocity derivatives include the
derivative of `i omega`, including at DC and for negative frequencies.
No damping, gyro coupling, source phase or nonsymmetric bearing coefficient
is removed, and no finite difference or modal truncation is introduced.

`RotorOperatingPoint::driven_resonance` and `survey_driven_resonances` pass
these actual responses and analytic tangents to the existing `fs-phs`
resonance owner. They retain its evaluation budgets, cancellation, isolated
peak/bracket requirements and grid-resolution limitations. Fixed signed
physical displacement/velocity receivers retain interference; their peaks
may differ from one another and from the free poles. Optional half-power
width and response Q belong to that observable and are not decay or material
loss estimates. Shared `exp(+i omega t)` phasors need no conjugation.
`RotorBearingModel::synchronous_response_with_derivative` follows
`omega = order * Omega` with the total shaft-speed derivative, including
the additional `i omega G` contribution. `RigidRotor` exposes
`unbalance_response_with_derivative` and `unbalance_resonance`; both retain
the original station phase/moment maps and `F = U Omega^2`, including
`F' = 2 U Omega`. Zero speed retains zero source and derivative. Peak
searches use positive shaft speeds and interpret returned frequency/width
as that same 1x speed. They compare steady operating points, not a spin-up
trajectory. The independent Jeffcott oracle in `rotor_unbalance_resonance`
distinguishes the source-dependent peak from both a fixed-force peak and
the natural frequency, and checks both analytic half-height crossings.
The original constant-coefficient linearization remains: a particular
solution does not establish decay of unstable homogeneous transients, and
varying bearings/preload/material state requires its own coefficient model.
G1/G3 tests in `rotor_resonance` cover analytic signed-frequency/source-phase
tangents, complete nonsymmetric/gyroscopic direction derivatives, distinct
displacement and velocity peaks, analytic linewidth, both whirl peaks, and
original binary detuning under coordinate rescaling. G4 cases retain
singular, budget, malformed-receiver and cancellation refusals.

`RotorOperatingPoint::driven_antiresonance` and
`RotorMaxwellOperatingPoint::driven_antiresonance` resolve coherent fixed-spin
receiver minima under the original fixed complex force/torque. They share the
native physical displacement/velocity projection and the `fs-phs` notch search.
Every retained mechanical coordinate and Maxwell stress contributes before the
signed receiver weights are applied; no observed response is inverted. A quiet
receiver can coexist with finite absorber motion, bearing stress and heat.
Response and tangent share one factor; recovery crossings skip unneeded tangent
solves. Original mechanical and constitutive residuals, source admission, query
limits and cancellation remain active. A caller-declared absolute recovery level
uses the squared original observation units, with both crossings required; its
width is not a free-mode or material Q. Three native tests check a tuned absorber
against an independent 2x2 inverse, a memory-bearing cancellation with finite
stress and one watt of branch heat, exact recovery-edge oracles, empty-memory
equivalence, invalid receivers, missing edges, cancellation and budgets.

`rotordynamics::periodic::periodic_floquet` forms the homogeneous period map
for 1–256 declared constant-coefficient portions of a physical cycle. Signed
spin and original full M/G are identical across the cycle; full C/K may change,
including anisotropic, cross-coupled and active coefficients. Coordinate order,
units and frame remain the caller's responsibility. Displacements and velocities
are continuous at coefficient switches. A changed stiffness may exchange energy
with its modulator; that exchange is not silently called heat.

A common invertible mass/time scaling retains all original coordinates, including
free modes. Native scaled Cholesky solves construct
`[0 I; -M^-1 K, -M^-1(C + Omega G)]`; independent double-double contractions of
the original mass and separate C/G terms gate each solved column. The existing
bounded `fs-phs` exponential advances each constant segment, and chronological
left multiplication forms `E_last ... E_first`. Native cancellable `fs-la`
complex QR returns all 2n multipliers. The report includes the original
`[q, q_dot]` monodromy, multiplier spectrum, spectral radius, numerical
`ln(radius)/period` growth rate and maximum original mass-equation residual.

There is no stable/unstable Boolean or rounding deadband around unit radius.
Any nonfinite or zero computed multiplier refuses as unrepresentable or
numerically unresolved, including loss of a fast-decaying direction while other
directions retain finite radius. No truncated state, artificial damping, pole
shift, symmetric bearing replacement or eigenvalue-error certificate is inferred.
This is the specified segmented model; using it to approximate smooth periodic
coefficients additionally requires coefficient/time refinement. Spin acceleration,
changing inertia, impulses at switches and nonlinear whirl are outside this API.

The existing 32-coordinate ceiling remains. Explicit limits admit segment count,
live dense bytes, whole-cycle multiply/add work and QR sweeps before dense
allocation. Cancellation reaches assembly, original-column residuals, matrix
products and native QR. Each shared exponential is indivisible but bounded to
64 states, 32 Taylor products and 48 squarings, with polls before and after;
no wall-time deadline is promised and no partial response is published.

Focused regressions include a closed-form parametric instability: two individually
conservative scalar segments with K/M of 1 and 4 and durations pi/2 and pi/4
have the period map diag(-1/2,-2). They also cover signed gyroscopic whirl,
full-inertia pole/exponential agreement, volume contraction, coordinate congruence,
cycle-origin invariance, free coordinates and unrepresentable/cancelled refusal.

`rotordynamics::periodic::{periodic_maxwell_floquet, prepare_maxwell_cycle}`
extends the homogeneous segmented cycle with EVERY original generalized Maxwell
stress. The physical state is `[q, q_dot, sigma_0, ...]`; supplied K is relaxed,
and branch stiffness must not also be added as an instantaneous spring. M, G,
signed spin and the original branch strain/stiffness/time remain fixed. Full C/K
may switch while displacement, velocity and stress remain continuous. All `2n+b`
modes survive, including equal-rate dark combinations and zero-strain memories.

One shared full-mass factor, common `[S q, T S v, T sigma/sqrt(k)]` coordinates,
cached memory columns, native segment exponentials, chronological composition
and native QR produce the full original monodromy and numerical spectrum.
Independent DD reconstruction checks original mass and memory-force/constitutive
rows. Full-state byte/work admission precedes dense allocation; cancellation
publishes no partial result. Prepared warm advances carry each returned physical
boundary into the next cached segment and audit original cycle increments.
Zero/unrepresentable multipliers refuse; near-unit rounding is not classified.
No heat/work ledger, eigenvalue enclosure, smooth-modulation error bound or
nonlinear bearing stability claim is implied by this motion path.

Six focused native tests cover exact full-state eigenmotions, independently
refined original-ODE RK4 modulation, parametric growth despite stable frozen
poles, dark/uncoupled memories, full-mass/signed-gyro coordinate covariance, all
66 original states beyond the mechanical-only ceiling, budgets and cancellation.

`rotordynamics::periodic::prepare_forced_cycle` extends the SAME preparation
owner with original generalized force/torque held over each declared segment.
It retains full M/G/C/K and the common mass/time coordinates. A normalized
constant source coordinate in the existing exponential yields each affine map;
it carries no physical mode or stored energy. Load amplitude therefore does not
set the exponential's scaling budget. The largest forced segment exponential
has 65 coordinates. Zero forcing uses the identical homogeneous exponential.

`PreparedRotorCycle` caches the full chronological affine map
`x_next = P x + b` and each segment's transition/source offset.
`advance` returns original physical q/v at every boundary, including the
initial state and complete period endpoint; the next segment starts from the
actually returned original state. Held force changes introduce no impulses.
Warm trajectories use matrix/vector operations without new exponentials or
factorizations. A rowwise DD check reconstructs the original increment equation
`(x_end-x_start) - (P-I)x_start - b = 0`. The unchanged state cannot hide a
small force behind its large magnitude; underresolved endpoints refuse.

`periodic_orbit` uses native partial-pivot LU on `(I-P)x=b` in the common
scaled coordinates. Its existing native condition estimator supplies
`||(I-P)^-1||_1 (1+||P||_1)` as a numerical resolvent-gain estimate with an
explicit caller limit. This is not an inverse-norm bound or an invariant
mixed-unit physical condition number. Singular maps refuse even when forcing
is compatible with a nonunique periodic family. An unstable cycle may have a
unique periodic particular orbit; no attraction or attained steady state is
inferred. The actual original initial state passes a rowwise force-scale DD
fixed-point check, and the independently propagated endpoint must close in the
declared common scaled full-state norm, including static/nodal zero coordinates.

Explicit byte/work budgets include retained segment maps and a complete boundary
trajectory/solve; the default dense-byte budget is 4 MiB. Cancellation publishes
no partial plan, trajectory or orbit. These are kinematics and numerical
residuals, not state-error certificates or full actuator-energy receipts.
Stiffness switches can exchange parametric work. Smooth forcing is a separate
refinement problem; held inputs are not silently reinterpreted as sinusoids.

Focused forced-cycle regressions cover analytic resonant square-wave response
across force scales, nonbinary static equilibrium, an unstable quarter-period
particular orbit and growing perturbations, full mutual inertia, coordinate
congruence and cycle origin, singularity, small-force resolution and cancellation.

`rotordynamics::periodic::power` adds physical work to original held-load,
continuous harmonic and finite Fourier cycles through
`prepare_forced_power_cycle`, `prepare_harmonic_power_cycle` and
`prepare_fourier_power_cycle`. Native maps and periodic-orbit solves still own
motion; the shared `fs-phs` stable Gramian kernel independently integrates
applied `f.v`, signed `v.C.v`, and nonconservative `v.K_skew.q`. Harmonic
source normalization, original-mass solves and common-origin phase preparation
are shared with native motion. Phase continues through coefficient and authored
amplitude switches; no held approximation or steady-state assumption is added.

The work integral observes DC and ALL simultaneous tones in one augmented
flow. This retains DC-tone, tone-tone and transient interference in every
quadratic work channel, including destructive interference in actual damping
loss. Separate per-tone powers cannot replace this integral. Cold workspace is
quadratic and integral work cubic in at most `2n+1+2*tones` coordinates (129
at the existing bounds), explicitly charged to the whole-query budget.
After integration, double-double contraction of the fixed source phases
preserves their complete cross blocks in one affine coordinate. Retained forms
and warm evaluation remain O(segments*n²), independent of tone count.

Polarized original endpoint quadratic forms independently measure kinetic and
symmetric-K storage increments. Every K switch, including the final-to-first
period seam, retains `q.(K_next-K_old).q/2` as parameter work at unchanged
state. The full-cycle balance uses first-segment storage on both sides;
average watts divide each physical work channel by the complete period.
This resolves actuator and stiffness-modulator consumption without adding
another physical trajectory solver.

Signed active-bearing and circulatory-stiffness work is never heat; exact skew G
does no work but remains in the motion equation. Native mass, residual, resolvent
and cancellation gates remain active, with whole-query budgets for all cached
forms and cold source/integral work. Each segment, accumulated absolute segment
defects and the independent whole-cycle energy residual pass the authored
energy allowance before publication. G1/G3 regressions target analytical damping
loss, DC/tone constructive interference, destructive damping interference,
signed gyroscopic forcing, global-phase subdivision, active C, circulatory K,
source/modulator exchange at the period seam, and resource/cancellation refusal.
No smooth-modulation, peak-power, physical heat routing, positive-storage,
attraction or certified state-error claim follows.

`rotordynamics::periodic::prepare_harmonic_cycle` drives those SAME physical
segments with continuous `Fc*cos(omega*t) + Fs*sin(omega*t)` forces/torques.
The positive integer harmonic order (at most 1024) fixes
`omega=2*pi*order/period`; phase is relative to the common cycle origin and
continues through every coefficient/load switch. The sine sign expresses
rotation direction independently of the signed shaft speed. For the existing
`exp(+i omega t)` force phasor, Fc is its real part and Fs minus its imaginary
part. Synchronous forcing requires the caller's period/order to match spin.

Two normalized source coordinates form a skew rotation inside the existing
bounded exponential (at most 66 coordinates). They carry no physical energy
and are excluded from the physical monodromy and Floquet spectrum. Both
original force columns pass the existing full-mass residual check. Source
amplitude does not set the exponential's scaling; a zero load uses the identical
homogeneous flow. No dynamic-stiffness inverse or held-sinusoid approximation
is introduced, including at a lossless resonance. The same cached trajectory
and unique-orbit solve retain singularity/resolvent and cancellation refusals.
A single harmonic and segmented coefficients do not imply smooth coefficient
modulation, spin acceleration or certified state errors; the companion power
API supplies the separate physical work integrals.

Six G1/G3 regression cases cover continuous resonant response across force
scales and subdivision, signed gyroscopic rotating loads, full mutual inertia
with free coordinates, secular lossless resonance, globally phased coefficient
switches, and homogeneous/admission equivalence.

`rotordynamics::periodic::prepare_fourier_cycle` composes steady physical
loading with up to 32 simultaneous harmonics of the same declared period.
Orders are strictly increasing in 1..=1024; every segment supplies its actual
DC and cosine/sine coefficients. All tones retain the common phase through
coefficient and amplitude switches. Empty or silent harmonic sets reuse the
held-source behavior, including the complete physical Floquet spectrum.

One shared mass factor and segment generator serve every component. The
existing normalized harmonic exponential integrates each tone separately,
and double-double accumulation combines their affine offsets. A single
physical transition is retained per segment, so warm trajectory/orbit work
and retained storage are independent of tone count. Cold exponential work is
linear in tone count with a fixed 66-coordinate maximum; the complete query
budget reserves all tones before allocation. Original full-mass force-column
residuals, cancellation and singular periodic-orbit refusals remain active.

Six G1/G3 regression groups cover analytic resonant multi-harmonic orbits
across 160 decades of force scale, full mutual-mass integration, globally
phased coefficient/amplitude switches, lossless secular growth, exact
held/silent-source parity, and whole-spectrum resource/shape admission.
This finite Fourier load does not infer an unmodeled force spectrum,
nonlinear bearing/mesh feedback, smooth modulation or spin acceleration;
waveform truncation remains separate from the companion finite-model work audit.

`harmonic_response` uses `exp(+i omega t)` and native cancellable complex LU
on `K - omega^2 M + i omega(C + Omega G)`. Row equilibration precedes the solve;
compensated inertial products retain the original binary coefficients' detuning
at narrow resonances, including mutual inertia and signed driving frequencies.
The accepted residual and signed mean powers use the original separate matrix
actions. Singular resonance refuses without a damping floor or pole shift.
Mean source, bearing-velocity, stiffness, gyroscopic and inertial powers are
reported separately. Cross-coupled/active work is not called heat. Equilibrium
residuals are measured numerical defects, not response-error certificates.
The harmonic solution is a sinusoidal particular solution: unstable homogeneous
transients need not decay, so it does not imply an attained steady response.

`RotorOperatingPoint::midpoint_stepper` now advances actual displacement/angle
and velocity states through the complete constant-speed rotor equation. Native
real partial-pivot LU factors `M + h/2 (C + Omega G) + h²/4 K` once; scaled mass
Cholesky also supports instantaneous full-inertia generalized impulse kicks.
Preparation costs O(n³); subsequent trajectory steps and kicks cost O(n²),
within the existing 32-coordinate bound. The step solves for a velocity
increment before updating the original state. The original impulse equation
is independently checked on the actual returned endpoint, with no mass lumping,
bearing symmetrization, damping floor or clipped impulse.

The supplied midpoint load drives the reusable trajectory state. Separate
signed input, bearing-velocity, stiffness and gyroscopic work accompanies each
step; active/cross-coupled bearing work is retained, and gyroscopic work is zero
in exact arithmetic. Symmetric conservative stiffness gives the usual midpoint
quadratic-energy identity, while phase accuracy still needs timestep refinement.
The physical rigid-rotor specialization preserves transverse angular momentum
`H = Jd theta_dot + Jp Omega [theta_y, -theta_x]` during unsupported constant-spin
whirl. The caller supplies load samples, coefficient validity and any event
time. Spin acceleration and nonlinear bearing evolution require their own
physical operators; the rigid specialization below supplies an explicit
prescribed-spin owner. Cancellation or failed balance publishes no partial state.

The transient tests exercise full-mass impulse-to-trajectory response, native
rigid-rotor angular momentum and energy, analytic whirl and damped ring-down
refinement, nonsymmetric active work, and singular/underresolved/cancelled refusal.

`rotordynamics::rigid::RigidRotor` supplies a concrete physical assembly. It
lowers the original full centroidal inertia by the supplied proper reference
rotation and admits axisymmetry under an explicit tight tolerance. Its four
stationary coordinates are lateral translation and tilt; finite bearings use
their actual axial stations through `A^T B A`, retaining force moments.
Signed-speed unbalance uses `A^T F`, and recovered bearing work keeps its sign.
This is a rigid small-tilt linearization, without flexible-shaft or
asymmetric-rotor Floquet claims. The steady lane retains constant speed;
the explicit ramp lane below also supplies the acceleration-dependent terms.

`RigidRotor::prepare_unbalance` retains a separate transient physical load
source for prescribed phase, signed speed and angular acceleration. Each
station supplies `U*(Omega²*e_r - alpha*e_theta)` and the original axial
force-to-bending-moment map. The tangential reaction remains present at zero
speed during startup and changes sign under braking. The original nominal
mass and inertia are unchanged. Phase is supplied independently from speed;
the caller integrates the actual speed history and samples loads at the
time-integration owner's evaluation times.

Two compensated complex source moments give O(stations) preparation and
fixed-size, allocation-free load queries. `force_vjp` returns analytical
phase/speed/angular-acceleration cotangents at fixed source geometry, sharing
the original load admission. Nonrepresentable assembled or phase-projected
moments, products and cancellation refuse without publishing a partial source.
This load source supplies neither a spin drive nor the additional physical
operators needed for a general variable-spin tilt or flexible-shaft model.
`tests/rigid_rotor_acceleration.rs` checks startup/braking forces and moments,
steady-phasor and phase-origin equivalence, all three load derivatives, native
four-coordinate translational response against an integrated acceleration
oracle with second-order refinement, and cancellation/representability refusal.

`RigidRotor::spin_ramp_step` advances those original four coordinates through a
prescribed constant-angular-acceleration interval. `RotorSpinState` retains
continuous time, phase and signed speed; phase is integrated quadratically and
the unbalance is sampled at its actual temporal midpoint, retaining tangential
forces during startup, braking and reversal. Prepared unbalance moments require
no station traversal during a step. A bounded four-coordinate native midpoint
factorization retains complete `M,C,K,G`, including anisotropic, nonsymmetric
and active bearing entries. The resulting trajectory is second order under
time refinement; constant acceleration zero reduces to the original steady
midpoint step without a second integration kernel.

The drive-torque direction is an explicit physical input. Differentiating
`H_perp = Jd theta_dot + Jp Omega [theta_y,-theta_x]` contributes both
`Omega G theta_dot` and `alpha G theta`. A drive following the actual tilted
shaft supplies the latter transverse torque and cancels it. A drive held
along nominal stationary +z leaves the full additional `alpha G theta`
operator. No torque direction is silently selected. Original bearing-stiffness
work and signed work against this acceleration coupling are reported separately;
their sum is the native effective-stiffness work. Acceleration coupling is
exchange with the prescribed spin boundary, not damping heat or elastic storage.

This owner prescribes spin rather than solving a finite-torque motor, drive
storage or flexible-shaft problem. It retains the original small-unbalance,
small-tilt and authored bearing-coefficient validity. A nonadvancing clock,
nonrepresentable kinematics/operator, native solve/balance refusal or cancellation
returns neither a partial lateral state nor a partial spin advance.
`tests/rigid_rotor_spin_ramp.rs` checks exact free angular-momentum precession
and its second-order refinement, shaft-following versus stationary torque,
analytic accelerating/braking/reversing unbalance translation, native steady
equivalence, independent bearing-point work and acceleration exchange, and
invalid input or final-poll cancellation with exact retry.

`RigidRotor::current_spin_motor` connects a balanced axisymmetric rotor's
original polar inertia to the existing finite-current electrical drive. The
shaft angle is its own physical coordinate, never one of the four lateral
translations/tilts. Shaft-axis actuation must be selected explicitly; nominal
stationary-axis actuation refuses. Nonnegative authored shaft drag and external
shaft torque enter the original scalar momentum equation. Native current,
terminal-voltage, back-EMF, speed-domain and electrical-source limits determine
the accepted current and actual spin acceleration, including stall and reversal.

The resulting native midpoint spin enters the original full lateral
`C + Omega G` operator. Native shaft displacement supplies continuous phase.
There is no prescribed acceleration, copied inertia parameter, eigensolve or
second integration kernel. Scalar and four-coordinate preparations remain
bounded; the smooth trajectory is second order under time refinement.

Spin kinetic storage, nonnegative shaft-drag work, signed lateral bearing work,
gyroscopic work, independent external mechanical work, electrical draw/return,
winding heat and bridge heat retain their original owners and meanings. Their
combined original energy balance must satisfy the authored power allowance
before any complete state is returned. Invalid lateral input or a final
cancellation also publishes no partial electrical debit or spin advance.

This owner retains the physical small-tilt truncation and authored bearing
coefficients. Residual unbalance is not admitted without its reciprocal spin
reaction; finite inductance, saturation and thermal changes are not inferred.
`tests/rigid_rotor_spin_drive.rs` exercises analytic power-limited acceleration
using the actual polar inertia, changing-spin free whirl and second-order
refinement, stall heating, reversal/regeneration, passive shaft drag, exact
native bearing-operator reuse, unsupported-axis/coordinate/domain refusals and
late cancellation with deterministic retry.

`ElasticEccentricRotor` supplies the finite-rotation reciprocal counterpart for
an actual eccentric rigid body on conservative elastic radial bearings. The
original `SpatialInertia` retains mass, COM offset from the authored shaft frame,
and the complete centroidal tensor. Native `MassProperties` are prepared once
at that COM; the offset is used only in material shaft-point geometry. Each
bearing force and its moment `arm_body cross R^-1 force_world` act together,
so lateral bearing loads accelerate or brake the actual spin. No prescribed
unbalance force, copied polar inertia or additional eccentric mass is inserted.

In the planar specialization, with shaft position `q`, COM offset `e e_r`,
`e_t = [-sin(phi), cos(phi)]` and COM polar inertia `J`, this reproduces
`m q_ddot + m e e_t alpha = F + m e Omega^2 e_r` and
`(J + m e^2) alpha + m e e_t dot q_ddot = shaft_torque`.
The parallel-axis term and both reciprocal acceleration couplings arise from
the original COM equations. Native point-motion queries retain tangential,
centripetal and full-inertia gyroscopic acceleration, including during reversal.

The prepared reference transform fixes bearing centres and stationary radial
axes; stations are measured from its authored shaft origin. Full finite pose
updates use the existing native configuration-force Verlet composition and
native free rotation. Bearing stiffness must be finite, exactly symmetric and
resolved positive semidefinite, and bearing damping must be zero. The original
matrix is retained without symmetrization or a support floor. Velocity damping,
active/fluid-film bearings, clearance, shaft flexibility and electrical source
limits are outside this elastic mechanical owner. The balanced electrical drive
above remains a separate owner; this API takes an explicit finite shaft torque
following the actual body +z axis and never prescribes speed.

Each step reports original native work, bearing storage and its polarized
increment, independently integrated endpoint-trapezoidal shaft work, and the
signed defect `delta(K + gravity potential + bearing potential) - shaft_work`.
Native bearing work is not counted a second time as source work. This is a
second-order method with a numerical energy defect, not exact conservation,
heat, or unconditional stability; users must resolve bearing and spin timescales.
State and load queries allocate no per-station scratch after bounded preparation;
cancellation or refusal publishes no partial candidate. The five regressions in
`tests/rigid_rotor_eccentric.rs` cover an independent coupled mass solve,
finite-torque reversal, planar ODE refinement with changing spin, original full
inertia and finite-rotation bearing virtual work, and refusal/late cancellation.

`ElasticEccentricRotor::current_spin_motor` binds the same original current,
terminal-voltage, back-EMF and electrical-source owner to the actual eccentric
full body. Shaft-axis actuation and scalar source coordinate zero must be
explicit. For a source impulse `j e_z` at fixed pose, the original tensor gives
`delta omega_z = j (I_body^-1)_zz`; its effective scalar source inertia is
`1/(I_body^-1)_zz`, not generally `I_zz`. The body's other angular momentum
components stay fixed during that source, while transverse angular velocity
and complete tensor kinetic coupling remain in the native impulse receipt.
This scalar inertia applies only to electrical source response and never
replaces the original body's trajectory inertia.

Two electrical half sources surround the original elastic native trajectory.
Each source reads actual incoming body-frame shaft speed, selects a current
through the original voltage/current/power gates, and applies its accepted
shaft impulse exactly once through `RigidBodyState::apply_momentum_impulse`.
Actual angular transfer, actual endpoint shaft rate and full-body kinetic work
must agree with the source within explicit allowances. Actual endpoint speed,
terminal voltage and source power also pass the original drive gates independently
of those transfer tolerances. Swallowed angular impulses refuse. The middle step
owns all physical pose motion, reciprocal bearing reaction, full-inertia free
rotation, gravity and separately supplied
external shaft torque. Scalar source displacement is solely a work quadrature;
it is not a second shaft angle or an extra elapsed physical interval.

The final balance reduces the original body's full endpoint kinetic increment,
gravity and bearing storage changes, original winding/bridge heat, source
return/draw, and external shaft work. Neither numerical mechanical error nor
source transfer error becomes electrical heat. Source draw and return remain
separate, including reversal, and every original electrical receipt is retained.
Explicit total-energy allowance and final cancellation admission precede return.
Preparation and a complete step use bounded scalar native solves and the
existing body/bearing owner, without a new nonlinear solver or inertia model.

The symmetric composition is second order for smooth held-current operation.
Active curtailment, controller transitions and continuous-time voltage/power
peaks require time refinement: electrical source-stage endpoint limits are not
a bound over the intervening nonlinear body trajectory. Finite inductance,
battery storage, thermal feedback and velocity-dependent bearings remain
outside this owner. Six regressions in `tests/rigid_rotor_eccentric_drive.rs`
cover original nonprincipal mobility and tensor work, analytic power-limited
startup, stall heat and regeneration, reciprocal eccentric-orbit refinement,
zero-current native equivalence, final cancellation and unresolved body impulse
refusal.

`DampedEccentricRotor` extends the same actual eccentric body with stationary
passive radial bearing damping. The original K matrices stay with the elastic
owner; each original C matrix passes the existing `PassiveBushingMatrix`
symmetry/positive-semidefinite admission and is retained without projection,
damping floors or diagonalization. Symmetric anisotropic cross terms, multiple
stations and finite body rotation retain their reciprocal COM moments. The
full source inertia remains the original native body tensor, including its
off-diagonal components; no housing body or additional rotor mass is introduced.

At a fixed actual pose, the source rates are
`v = [v_COM_world, omega_body]`, with `M = diag(m I3, I_body)`.
Each stationary radial axis has the conjugate material-point velocity row
`J = [e_world, arm_body cross R^-1 e_world]`. All bearings are assembled
simultaneously as `D = sum J^T C J`; the original full-mass `RotorMidpoint`
owner solves `(M + h D/2) delta_v = -h D v`. One original native momentum
impulse applies the combined force and central moment. This is a source at
fixed pose; native free rotation owns gyroscopic transport. The source's
auxiliary displacement is a work quadrature and never a second physical pose.

Each accepted source checks actual native endpoint rates and momentum transfer
against that original constitutive law with explicit velocity/impulse budgets.
Per-bearing heat is independently `h u_mid^T C u_mid` at the actual returned
body midpoint rates, evaluated with compensated original-coefficient products.
Both actual kinetic change plus heat and original impulse work plus heat must
pass an explicit energy allowance. Heat is not obtained by relabeling kinetic
loss or a numerical residual. Nonfinite or unresolved negative quadratic work
refuses; no clipping turns an unverified response into passive dissipation.

The source enforces a conservative original-inertia decay-rate trace bound
`trace(M^-1 D)` with a numerical roundoff margin. The caller's maximum
source-duration times that bound must lie in `(0,2]`. Implicit midpoint is not
L-stable: unresolved stiff modes would otherwise change sign while apparently
losing energy. This gate prevents that modal reversal and requests refinement;
it is neither an interval spectral certificate nor a full trajectory accuracy
estimate. Null damping directions remain free, and zero damping or duration
skips the source solve.

A complete step composes a damping half source, the original elastic
configuration-force Verlet trajectory, then a damping half source at its
actual endpoint pose. No velocity force is passed to configuration-only Verlet.
For smooth fixed coefficients and resolved elastic/free rotation this symmetric
composition is second order. The final original-body balance is
`delta(K + gravity potential + bearing potential) + bearing_heat - shaft_work`.
Physical bearing heat, external shaft work and signed numerical mechanical
error remain separate, with a whole-step energy gate before publication.
Instantaneous load and native point-acceleration queries also include damping
and its reciprocal shaft torque. Refusal/cancellation returns no partial body
candidate, and the public fixed-pose source can be composed by a coupled owner.

This owner admits stationary bilateral supports and a held finite body-axis
shaft torque. Active/nonsymmetric bearings, moving supports, clearance, inferred
speed/temperature coefficient laws, fluid-film states and a thermal destination
require their original constitutive owners. Electrical and passive source laws
retain separate receipts; the explicit current-drive composition below combines
them. Calling only the elastic component does not add damping. Five regressions in
`tests/rigid_rotor_damping.rs` specify scalar reciprocal spin braking and analytic
midpoint heat, full-tensor/multiple-bearing dense equations and virtual power,
independent planar ODE and heat refinement, passive/modal refusal, zero-damping
native equivalence and final cancellation with deterministic retry.

`DampedEccentricRotor::current_spin_motor` binds the original current-driven
eccentric owner to that same passive body, making electrical drive and bearing
damping one admitted physical trajectory. Explicit shaft-axis actuation remains
required. `DampedEccentricDriveControl::drive` retains the original electrical
limits, transfer checks and middle-trajectory energy gate; its duration is the
sole physical interval. The original damping controls retain their material,
source-transfer and modal-resolution limits, and their total-energy allowance
now gates the complete coupled trajectory.

The composition is `D(h/2), E(h/2), elastic(h), E(h/2), D(h/2)`.
The existing electrical owner supplies both E sources; the existing passive
owner supplies both D sources. Every source reads actual incoming body momenta.
Both owners borrow the same original elastic geometry and full tensor; the
electrical effective source inertia never replaces trajectory inertia. Only
the original elastic stage advances pose, while the source intervals integrate
their respective physical work and heat without adding elapsed physical time.
No extra body, duplicate force integration or new controller solve is introduced.

The final passive source can increase body-z speed by transferring COM motion
into spin even while reducing total kinetic energy. Therefore the middle
electrical receipt cannot establish the final motor operating point. The
composition checks actual initial and final shaft speed, terminal voltage and
signed source power with the adjacent actually accepted first/last current.
Both draw and regenerative-return limits remain authoritative. Endpoint failure
returns no complete candidate or electrical debit; it does not silently
reselect a current after its work was integrated. These original-law endpoint
diagnostics add no extra source energy.

The outer ledger recomputes original initial-to-final full-body kinetic change
and retains original gravity/elastic storage increments. It reduces the
original half-source winding heat, bridge heat, draw and return separately,
adds independently measured bearing heat, and subtracts only the actual
external shaft-work receipt. Electrical shaft work is internal exchange and
is not counted again as an external source. No signed residual is converted
to physical loss or a controller efficiency. The complete energy gate and
final cancellation check precede publication.

This symmetric source composition is second order for smooth fixed
coefficients and resolved mechanics with held uncurtailed current. Active
curtailment, controller switches and nonlinear interstage electrical peaks
require refinement; endpoint checks are not continuous-time peak certificates.
The underlying stationary-support, current-regulated zero-inductance and
fixed-parameter material restrictions remain explicit. Battery storage,
inductance, thermal feedback and moving-support work require their own owners.

Four focused regressions in `tests/rigid_rotor_damped_drive.rs` cover an
independent powered/damped planar acceleration and ODE/energy-integral
refinement, successful source curtailment and regeneration with all three
physical losses separate, final damping that spins up the rotor beyond an
otherwise accepted middle drive's speed/voltage/draw/return limits, both
zero-source reductions, and late cancellation/whole-energy refusal.

`DampedEccentricRotor::voltage_spin_motor` supplies finite-inductance voltage
actuation on that same actual eccentric body. `EccentricRotorWindingState`
carries the sole native body and signed winding current across both electrical
half sources and into the next complete step. Magnetic storage derives from
`L I^2 / 2`; current is never reconstructed from a resistive equilibrium,
clipped after limiting, or reset between mechanical stages. Explicit shaft-axis
actuation and strictly positive native winding inductance remain required.

The parent's original full-tensor shaft-mobility source is reused unchanged.
The existing `RotorVoltageMotor` solves its coupled midpoint winding and
mechanical equations and selects the feasible held terminal voltage under
original voltage, current and bidirectional source-power limits. It prepares
one cached coupling per complete step and reuses it for both half sources.
The trajectory is `D(h/2), E_L(h/2), elastic(h), E_L(h/2), D(h/2)`;
purely mechanical stages leave the actual current unchanged, and only the
original full-body elastic stage advances pose. No extra rotor, replacement
inertia, winding solver or iteration is introduced.

Each accepted midpoint torque is applied once through the original body
momentum-impulse owner. Actual angular transfer, actual body-z rate and
complete tensor kinetic work must agree with the source within explicit
allowances. The winding equation is then independently re-evaluated using
actual returned body midpoint speed and actual current endpoints, so a scalar
prediction tolerance cannot waive the original electrical voltage balance.

Actual terminal voltage is the source's accepted V, and terminal power is
`V I`. The zero-inductance equation `V = Kt omega + R I` is not a voltage
compliance condition here: `L dI/dt` supplies the remaining voltage. A shorted
spinning rotor can therefore develop braking current and magnetic storage with
zero terminal voltage and zero source draw. Outer speed/current/terminal/source
admission uses the actual complete endpoints after bearing damping, while the
current and accepted V retain their original meanings. Source curtailment
that leaves no current-feasible operating point refuses atomically rather
than erasing current or its stored magnetic energy.

The complete ledger recomputes original full-body endpoint kinetic change and
the polarized endpoint magnetic increment `L I_mid_total (I_after-I_before)`.
It retains original gravity and elastic storage, bearing heat, winding
midpoint Joule quadrature, bridge heat, electrical return/draw and external
shaft work separately. Electrical shaft work is internal transfer and is not
counted twice. Native `R I_mid^2 dt` remains midpoint heat quadrature, not an
exact RMS-current integral. No balance residual becomes physical heat.

For smooth fixed coefficients and uncurtailed held voltage, the symmetric
composition is second order with resolved winding, bearing and free-rotation
timescales. Midpoint is not L-stable; stiff transients and active controller
limits require refinement. Endpoint admission does not certify continuous
electrical peaks. Battery storage, bridge switching, thermal feedback,
saturation and changing material coefficients remain separate owners.

Six regressions in `tests/rigid_rotor_voltage_drive.rs` cover independent
short-circuit electromechanical equations and current continuation across
complete steps, full-tensor source mobility/work, finite-current startup and
powered damped planar ODE/energy refinement, source curtailment/regeneration
with back-EMF above terminal voltage, final damping speed-domain refusal and
late cancellation, and a source-cutback/current incompatibility that must
preserve existing magnetic state for deterministic retry.

Focused regressions use independent Jeffcott amplitude/phase, gyroscopic circular
mode splitting, dense coordinate congruence, cross-coupled/active work and
original force residuals. A manufactured nonsymmetric rotor has an exact
sign-changing marginal speed at 0.5 rad/s; a closed-form isotropic gyroscopic
rotor has two exact 1x and two exact 2x synchronous crossings. Those oracles
check bracket width, count direction, gyroscopic-sign invariance, bounded work,
cancellation and malformed-grid refusal without assigning modal labels. They do
not close the broader certified Campbell/Floquet/whirl bead.

## Spatial gravity on native rigid-body dynamics

`gravity::PointGravityField` exposes specific potential, acceleration and the
analytic world-coordinate acceleration gradient in one conservative sample.
The uniform and central implementations retain the same original field, frame,
domain and potential datum. Central derivatives use radial factors without
cubing radius or taking numerical differences; the point path allocates nothing.
Sample admission requires finite values and a symmetric gradient. It cannot
prove that a custom implementation differentiates one stationary potential.
`tests/gravity_field.rs` checks spatial/potential finite differences, zero vacuum
trace, rotations, the native rigid mass owner, frame/domain refusals and extreme
radius arithmetic. The derivative API refuses unrepresentable gradients without
changing the legacy acceleration-only admission domain.

`gravity::SuperposedGravity` borrows a caller-bounded list of stationary sources
in one declared inertial frame. It implements the existing point and/or rigid
field traits, including heterogeneous trait-object lists, so native rigid,
articulated and deforming-solid consumers retain their original mass, pose,
quadrature and integration owners. Each component's potential, force, body
torque and analytic acceleration gradient are accumulated in deterministic
source order with double-double arithmetic. Sampling allocates nothing and
performs O(source count) work; existing consumers poll cancellation around these
bounded queries. Nested lists retain their own child budgets. Empty lists are
zero fields but still check frame, position and native mass admission.

Direct potential increments are combined before rounding the total. A custom
component without this optional capability is resampled at its two endpoints
for its own fallback difference, without discarding accurate work from the
other components. Callbacks must remain repeatable and side-effect-free;
fallback accuracy remains limited by the component's endpoint potentials.
Component refusals propagate without a partial load or candidate state. These
are prescribed stationary sources, not a self-consistent recoiling system or a
path-clearance guarantee. `tests/gravity_superposition.rs` checks potential/
force/tangent consistency, compensated balanced forces, mixed direct/fallback
work, native ballistics and cancellation, tidal torque, and frame/domain/budget
admission. `fs-couple/tests/solid_gravity.rs` exercises combined source forces,
Newton tangents and datum-independent work on original solid quadrature.

The optional `PointGravityField::specific_potential_increment` and
`RigidGravityField::potential_energy_increment` evaluate the same stationary
potential's endpoint change directly. Uniform fields cancel the reference
datum before compensated force/displacement arithmetic. Central fields use a
scaled rationalized inverse-radius difference, including transverse motion
too small to change a rounded radius. Custom fields default to `None`; owners
then subtract their already retained endpoint samples without another callback.
These scalar differences do not replace endpoint force/domain admission.
Rigid and point-source quadrupoles retain their correction's endpoint change
separately from the monopole, so small tidal work survives the total potential's
rounding. They reuse the original quadrupole correction and expansion domain.
Native and articulated gravity receipts retain the direct potential change for
work accounting while keeping the original absolute diagnostic energies.
The field, force stages and integrator are unchanged, and finite-step work
residuals remain measured numerical defects rather than dissipation.

Central and quadrupole sources also expose `specific_potential_increment_to`
between two explicit source snapshots and two sampled material points. Both
source and sample may move. Central snapshots retain G M and frame; quadrupole
snapshots retain the original mass, full inertia, G, extent and frame. The
shared rationalized radial increment includes source recoil before relative
coordinate rounding, and each original quadrupole correction retains source
orientation work separately from the monopole. Both endpoint domains remain
active; source reaction wrenches and path clearance belong to their trajectory
owner. Stationary increments delegate to this same native capability. New
field regressions cover moved-source and transverse work, common translation
with source rotation, reversal, isotropic reduction and changed-source refusal.

`gravity::CentralGravity` owns frame-checked Newtonian acceleration and specific
potential, with caller-supplied source position and G M. Rocket environments
reuse this owner. `UniformGravity` binds the existing acceleration to an explicit
potential datum. Neither infers Earth, a radius, softening or a rotating frame.

`gravity::UniformSphereGravity` adds an exact finite, uniform-density spherical
source with explicit center, radius and G M. Outside the source it reuses the
original central field by the shell theorem. Inside, acceleration is linear in
distance from the center and potential is quadratic; the center has finite
potential, zero acceleration and a finite isotropic gravity gradient. The point
sample's analytic acceleration gradient supplies the same potential's Hessian
to original deforming-solid quadrature and Newton tangents. At the declared
density discontinuity the force and potential are continuous, while the tangent
is explicitly the exterior one-sided value; no unique surface Hessian is claimed.

Direct potential increments retain interior and surface-crossing work without
subtracting two large absolute binding potentials. Existing superposition and
native adaptive rigid motion consume the same field interfaces. Queries perform
constant work without allocation. The rigid adapter samples the original target
mass at its COM and supplies no finite-target gravity-gradient torque; distributed
targets use their existing material quadrature. The spherical source is prescribed
and stationary, with no source recoil, deformation or inferred material profile.
`tests/gravity_sphere.rs` checks shell-theorem recovery, interior force/potential/
tangent consistency and Poisson trace, surface behavior, retained work and native
harmonic motion through the source center.

`gravity::RadialGravity` generalizes the prescribed source to concentric layers
of uniform density, each declared by its outer radius and shell mass. Positive,
increasing radii and nonnegative masses admit dense cores, lighter envelopes and
vacuum cavities; all-zero masses define a zero source. The explicit layer-count
budget is checked before allocation. Construction prepares compensated enclosed
gravitational parameters and exterior-shell potentials in O(layer count) work
and storage. Point samples locate their layer by binary search, allocate nothing
and evaluate the exact shell-theorem force, potential and analytic tangent.
Thin-shell arithmetic factors radius differences before forming ratios.

The field is nonsingular at the center, zero in force inside a hollow shell,
and reduces to the original central source outside the declared support. Density
interfaces preserve force and potential; tangents select the exterior-side local
density. Scalar acceleration/potential queries do not require a representable
tangent. Unrepresentable layer/total G M or central potential refuses construction.

Direct work uses the original enclosed mass and positive shell-volume terms,
splitting radial motion at crossed interfaces. Compensated endpoint displacement
retains small transverse work even when rounded radii and absolute potentials
are identical. Point work traverses only crossed layers, bounded by the admitted
profile; rigid work multiplies the same increment by the original target mass.
The COM-target, prescribed-source and continuum-quadrature semantics match the
uniform sphere. `tests/gravity_radial.rs` checks uniform-profile recovery, hollow
and dense-core analytical fields, Poisson/interface tangents, independent
potential derivatives, dominant outer-shell potential cancellation, cross-layer
work and reversal, thin layers, zero sources and input/frame/budget admission.

`gravity::QuadrupoleSourceGravity` exposes a native extended source's original
mass, full inertia and COM pose through the same potential/force correction and
its analytic spatial derivative. A declared inertia-compatible source radius
and radius/distance ratio admit far-field point and continuum samples. Isotropic
inertia reproduces `CentralGravity` exactly. Sampling a source snapshot supplies
no recoil or spin update; reciprocal consumers apply the equal/opposite force
and source COM torque from that same potential. No interior-gravity or
higher-multipole accuracy is claimed.

`gravity::step_in_field` uses the native `step_verlet_with` endpoint-load map
for conservative force/torque half-kicks and samples the actual
endpoint pose. It returns the native state/diagnostic receipt with the same
field's potential. Translation is velocity Verlet; no separate mass, attitude
or momentum owner is introduced. Uniform-gravity ballistics remain exact up to
roundoff; resolved Kepler motion has second-order phase/energy error. The
stationary conservative map is symplectic without an additional nonconservative
wrench; this is not exact energy conservation or an unconditional stability bound.
Cancellation is polled before sampling, before/after drift and before publication;
failures return no partially accepted state.

The explicitly selected quadrupole model uses the full existing COM inertia.
Its potential `-mu*m/r + mu*(3*n.I.n-tr(I))/(2*r^3)`, torque and translational
force are derivatives of one energy. The orbital force correction is retained
alongside gravity-gradient torque, preserving total angular momentum about the
source. No per-sample tensor decomposition or heap allocation is required.
The caller supplies a COM bounding radius and maximum radius/distance ratio;
necessary mass-radius consistency and endpoint limits are checked. This is a
second-mass-moment approximation, not a certified truncation error, resolved
finite-body gravity, self-gravity, source recoil, J2 or relativistic model.

`tests/gravity.rs` checks independent potential derivatives, full rotated inertia,
two-mass multipole refinement, exact uniform recovery, eccentric Kepler refinement
and long-time energy behavior, tidal angular-momentum balance, singular/frame/
extent refusal and late cancellation (G1/G2/G3). Endpoint samples do not certify
that a trajectory avoids the source or the model boundary between samples;
the caller must resolve close approaches and provide required event admission.

### Closed mutual-gravity systems with source recoil

`gravity::mutual::MutualGravity` advances slices of the original `RigidBodyState`
and `MassProperties` under direct Newtonian monopole interactions. Every unordered
pair supplies one potential `-G*m_i*m_j/r` and equal/opposite world COM forces.
Both source masses move. The caller supplies G in SI, one inertial frame, a
nonnegative separation cutoff and an explicit pair-count budget; no Earth value,
softening or far-field solver is inferred. Work is O(N²), storage is O(N), and
compensated force/diagnostic reductions follow the supplied deterministic order.
Binary-scaled pair coefficients avoid intermediate overflow/underflow for widely
separated mass scales; pair admission is unchanged when its bodies are exchanged.

All first half-kicks precede the original native force-free rigid drift; all
endpoint forces are then sampled together before the final half-kicks. Native
full-inertia spin is retained without a monopole torque. The receipt accounts for
translational/rotational kinetic energy, pair potential counted once, total linear
momentum and orbital-plus-spin world angular momentum. Smooth fixed-step motion
is second order and symplectic; energy is not exact, and bounded long-time behavior
requires resolved bounded trajectories. No momentum or energy projection is used.

`potential_energy_change_j` retains the direct sum of physical pair-potential
increments separately from those absolute diagnostics. The shared rationalized
inverse-radius owner uses both bodies' actual COM displacement before rounding
the relative separation, retaining small transverse work and source recoil.
The original binary-scaled pair coefficient remains authoritative. The negative
of this increment is gravity work, without implying exact finite-step energy closure.
This accumulation shares the existing drift-domain traversal and adds no
per-pair storage. Empty and single-body systems have zero potential change.

Pair work is admitted before allocation. Cancellation is polled per body/pair and
at publication; refused allocation, arithmetic, frame, budget or separation leaves
the immutable input system untouched and returns no candidate. Both endpoints and
the complete straight relative drift segments must clear the separation cutoff,
preventing a split drift from tunneling through a source. This is not a validated
enclosure of a true curved encounter: callers must resolve near-field motion.
The shared monopole/quadrupole drift guard retains both cross-product residuals
before exponent rescaling. Nearly collinear encounters therefore do not acquire
a fictitious clearance from independently rounded products; exact diagonal
crossings remain excluded. The regression includes all coordinate permutations,
large-distance scaling, excluded encounters and admitted near misses.
`tests/mutual_gravity.rs` covers analytic finite-source recoil, independent force/
potential derivatives, the two-body Kepler period and second-order apocenter error,
20-orbit energy/momentum behavior, frame/boost/permutation covariance, reversal,
native rotated-inertia spin, hidden interior source crossing and late cancellation
(G1/G2/G3). This monopole lane does not include finite-body effects. Self-energy
within a body, collisions, deformation, external loads, variable mass and FMM
remain outside this closed-system map. Adaptive interval advancement composes
this same step through the explicitly selected controller described below.

### Recoiling extended bodies and tidal spin/orbital exchange

The explicit `MutualGravity::quadrupole` extension borrows one declared COM
bounding radius per original body and reuses its full native inertia. For
`d = x_j-x_i`, `n = d/r` and `Q = m_j I_i + m_i I_j` in world axes, one pair owns
`U = -G*m_i*m_j/r + G*(3*n.Q.n-tr(Q))/(2*r^3)`. Both COM forces and both individual
gravity-gradient torques differentiate this potential. Each shape correction
reacts on its source; the same pair force is applied with opposite signs.
The noncentral orbital moment balances the two spin torques. No source is fixed
and no after-the-fact momentum correction is used.

The second-moment force, potential and torque arithmetic is shared with the
original fixed-source quadrupole owner. Native `begin_verlet`/`finish` stages
advance all bodies with simultaneous endpoint loads and the original full-inertia
free rotation. The existing system receipt carries both kinetic energies, pair
potential, linear momentum and total orbital-plus-spin angular momentum. The
resolved fixed-step map is second order and symmetric; energy has a timestep
defect and is not claimed to be exact.

The same direct potential-change receipt keeps each body's quadrupole correction
change separate from the monopole binding energy, using original admitted
inertia traces and the original correction owner. Small tidal spin/orbital work
therefore survives when the two rounded total potentials are identical. Focused
regressions independently check both rotating shapes' analytic energy changes,
small transverse pair motion, finite source displacement, translation/body-order
invariance and reversal without changing force or integration stages.

The sum of the two bounding radii divided by separation must not exceed an
explicit ratio strictly below one. Radius values, necessary inertia/radius
consistency, the original separation cutoff and pair budget are admitted.
Both endpoints and every closest straight split-drift separation meet the same
combined-radius bound. Overlapping declared bounds cannot pass. These are caller
shape declarations and an approximate far-field domain, not a certified support
or truncation-error bound or an enclosure of a true curved encounter.

Work remains O(N^2), scratch O(N), and pair evaluation performs no tensor
decomposition. Original compensated reductions, allocation refusals, per-body/
pair cancellation and all-or-nothing publication are retained. Derived
coefficients that cannot be represented refuse rather than discard a correction.
Quadrupole-quadrupole products, higher moments, resolved near-field mass geometry,
deformation and contact remain outside this opt-in model.

`tests/mutual_quadrupole.rs` adds independent translation/rotation derivatives,
full rotated inertia, body permutation, closed-system tidal momentum exchange
and reversal, a resolved 36-pair mass-cloud potential oracle with fourth-moment
refinement, combined-radius/interior-drift refusal, budgets and late cancellation.

### Exact mutual gravity of finite spherical mass distributions

`MutualGravity::uniform_spheres` binds a positive COM radius to each original
body. Every native principal inertia must agree with `2*m*R^2/5` within a
128-epsilon relative roundoff allowance; source density and rigid inertia thus
describe the same uniform sphere. Original masses, poses, momenta, frame, G and
pair budget remain authoritative. No replacement mass or rotational state is
introduced.

`MutualSphericalGravity::load` integrates the Newtonian interaction of both
uniform densities. Exterior pairs recover the original point-source owner;
full containment has a quadratic separation potential; partial overlap uses an
analytic spherical-cap integral. Force differentiates the same potential, with
one equal/opposite pair contribution and no COM torque. Coincident centers have
finite interaction energy and zero force. The pair law uses bounded dimensionless
overlap factors and scaled coefficients, rather than radius cubes or runtime
volume quadrature. Compensated sums retain deterministic input order.

`MutualSphericalGravity::step` uses the same native half-kick/free-drift/half-kick
stages as point mutual gravity. Both finite sources recoil under simultaneous
pair forces, while each original isotropic spin follows the native rotational
map with zero gravity torque. The resolved fixed-step map is second order and
symmetric, with a finite timestep energy defect; linear and orbital-plus-spin
angular momenta are conserved to roundoff. Direct pair-potential increments
retain both source displacements, including transverse work below the rounding
of either absolute binding energy. No extra pair traversal, per-pair storage,
replacement integrator, or energy/momentum projection is introduced.

`MutualSphericalGravity::advance` applies the original bounded step-doubling
driver to that fixed-step map. Accepted states and accumulated work come from
the two native half steps. Body-pair and trial budgets both remain binding;
invalid inputs, arithmetic refusal, exhausted budgets and cancellation publish
no partial state or interval. Local error control does not certify global error
or the continuous path, and variable-step advancement is not claimed symplectic
or time reversible.

In this finite-source model, a zero cutoff admits coincidence. A positive cutoff
retains the original strict exclusion at endpoints and throughout each straight
relative drift. Zero cutoff skips that unnecessary point-singularity guard and
allows center crossing. Source/radius/inertia/count/budget errors and cancellation
return no partial load. Overlapping density distributions are
mathematically admitted by gravity; contact and deformation remain with their
separate physical owners. The pair diagnostic excludes each body's constant
internal self-energy. `tests/mutual_spheres.rs` compares against independent
shell-volume quadrature, analytic containment/exterior/equal-radius limits,
potential derivatives, coordinate/body-order covariance and native admission.
`tests/mutual_spheres_dynamics.rs` checks center crossing, analytic containment
orbits and second-order refinement, direct small work and regime transitions,
native spin, reaction/momentum balance, reversibility, adaptive accuracy and
bounded atomic refusals.

### Adaptive gravitational encounters on native states

`MutualGravity::advance`, `MutualQuadrupoleGravity::advance` and
`MutualSphericalGravity::advance` advance a requested interval using the existing
force, recoil, full-inertia drift and potential-work owners. Each trial compares
one full step with two half-steps. The accepted state
is the actual two-half-step result; it is never Richardson-extrapolated or
projected onto an energy/momentum constraint. The second-order local estimate is
one third of their difference. Position, velocity, relative attitude and world
angular velocity each have an explicit absolute tolerance in SI units; relative
scales use local changes, so an arbitrary world origin or common velocity does
not inflate the requested error allowance. Every body participates in the norm.

`GravityAdvanceOptions` bounds the initial, minimum and maximum trial durations
and the total number of trials. Accepted and rejected trials share that budget;
each trial performs at most three original native steps. The result exposes the
counts, largest accepted error ratio and a suggested next duration. Only accepted
half-step potential changes enter a compensated interval work sum. Endpoint
diagnostics retain the native original mass, pose, momenta and potential owners.
Cancellation, invalid data, unrepresentable arithmetic, exhausted budgets or an
unresolved minimum step return no candidate interval and leave inputs untouched.
Initial domain admission is mandatory; only an explicit finite-body expansion
domain or pair-separation cutoff refusal during a trial may be refined. The
typed `GravityError::SeparationCutoff` retains the rejected separation and its
declared strict lower bound. Other native/source errors retain their identity
and are not silently retried. The original drift guard remains active at every
trial size; retrying does not authorize a trajectory inside the excluded region.

The estimate is local and asymptotic, not a certified global error bound or a
continuous-path enclosure. Native cutoff and expansion-domain checks remain
binding; the controller is not a collision/event solver. State-dependent step
selection is not claimed symplectic or time-reversible, even though each native
fixed-step map retains its original structure. No discrete adjoint or bitwise
equivalence between differently partitioned interval requests is claimed.
`tests/adaptive_gravity.rs` exercises physical orbit accuracy/refinement, native
rotational dynamics, momentum/work accounting and bounded atomic refusals.
`tests/adaptive_gravity_encounters.rs` resolves a circular orbit whose coarse
split drift falsely crosses the cutoff, for both monopole and spherical
quadrupole models, and refuses actual cutoff entry and inadmissible initial data.

`gravity::advance_in_field` applies the same controller to the original
`step_in_field_with_loads` map for uniform, central, quadrupole, superposed and
custom stationary fields. The single state is retained in a fixed-size array,
so the adapter and controller allocate no heap storage for this path. Held
external world forces/body torques keep their original substep semantics;
their work is excluded from conservative potential changes. Source callbacks
must be repeatable and side-effect-free. Invalid options are refused before
source evaluation, and original frame/domain/arithmetic errors propagate.
Initial and endpoint samples do not certify the intervening trajectory.
`tests/adaptive_gravity_field.rs` checks analytic combined-source ballistics,
small direct work under a large potential datum, eccentric-orbit convergence,
deterministic replay and atomic late source/cancellation refusals.

## Configuration-dependent loads on the original rigid integrator

`RigidBodyIntegrator::step_verlet_with` composes two endpoint force/torque
half-kicks around the original momentum-preserving full-inertia free map.
The first load belongs to the supplied state; the callback evaluates the
actual drifted endpoint and returns its original load receipt with its
wrench. Gravity acts once through the native owner. Each body torque uses
its own endpoint axes, so its world impulse is rotated at that endpoint.

The staged owner retains the linear half-kick in double-double precision until
both the position contraction and endpoint kick finish. Representable resultant
impulses survive half-kicks below the stored momentum's rounding threshold;
cancelling endpoint forces also retain the original momentum. Position
is rounded after adding the original position, preserving a small displacement
when large drift and origin offsets cancel. Callback views and accepted states
keep their original `f64` representation; no extra substeps or allocation occur.

Angular half-kicks retain the same two-component momentum through the native
`fs-time::lie::compensated` Hamiltonian drift and the endpoint kick. Retained
momentum drives each rotation phase, quaternion composition and paired inverse
momentum transport; carrying a tail after a rounded drift is insufficient.
The cached principal frame supports full SPD inertia without a new eigensolve.
`step_symplectic` uses this same staged owner for its held wrench, while
`step_kick_drift` retains its full kick through the same angular drift and
linear position contraction. Half impulses include any subnormal-duration
remainder before force contraction, preserving a representable complete impulse.
Equal endpoint forces use one finite full impulse for accepted linear momentum,
including minimum subnormal impulses whose individual halves cannot be stored;
overflowing full products retain the original staged evaluation.
Translation in `step` and `step_cayley` also uses the shared retained half-kick
and position contraction alongside their existing implicit angular midpoint.
The default `step` also retains the angular midpoint through its exponential
shortcut for exactly spherical or single-principal-axis motion. This uses the
same common free owner where the gyroscopic term vanishes; the Cayley chart
and general asymmetric nonlinear solve retain their existing equations.
`tests/angular_momentum_stages.rs` supplies independent large-phase, full-inertia,
endpoint-impulse, reversal, full-kick quotient, subnormal-duration and cancelling
linear drift regressions across the public rigid steppers. These are numerical
precision improvements within the existing splitting and accepted state types;
they add no energy projection, extra time subdivisions or general-inertia solve.

For smooth configuration-only loads this is a symmetric second-order
splitting. Spin plus orbital momentum changes only through the
two physical endpoint impulses, up to floating-point error. The free map
retains its existing exact spherical/axisymmetric cases and second-order
asymmetric-body energy defect. No energy projection or angular-velocity
rescaling is introduced. Timestep refinement remains required; this does
not establish second order for velocity-dependent loads, discontinuous
events or arbitrary adaptive steps. Callbacks stage external histories
until the composing owner accepts the complete result.

Native full-tensor tests check inertial force/torque impulse balance,
time reversal under configuration-dependent conservative loads, and late
source/wrench refusal without publication. The hydrostatic material owner
uses this same API for finite-rotation righting motion.

## Directional friction in configuration-admitted rigid steps

`frictional_impact::step::elliptic` applies the simultaneous directional
contact owner inside the original force-kick/contact/torque-free-drift engine.
Scalar and directional paths share the same initial/endpoint gap checks,
complete-path geometry admission, original pose integration, momentum balances
and work accounting. Normal/tangent and shared-body coupling remain intact.
Force-splitting work and drift defects remain numerical terms, not friction heat.

`step_rigid_elliptic_contacts_with` samples each active source from the actual
post-kick, precontact velocity. Private source candidates survive through the
native impact, drift, geometry callback and final whole-step audits; bodies and
histories publish together only after all succeed. Inactive histories remain
unchanged. Staging requires no clone bound on accepted source histories.

This retains the original first-order held-force splitting and its timestep
refinement requirement. The caller supplies the normal law, held friction role,
original collision/path authority and any objective texture transport.
Instantaneous loss does not define finite slip history, pressure, wear or heat
allocation. Nonnegative principal coefficients include exact free slip along
one principal tangent axis. Every response remains subject to native physical
mobility and determinacy admission.

## Directional friction through original articulated contact

`articulated_contacts::friction::elliptic` routes oriented principal
coefficient pairs through the original point Jacobians, articulated-body
inverse-mass actions and complete native endpoint/momentum/work audit.
Normal/tangent, joint/base and shared-contact coupling remain simultaneous.
No generalized dense mass matrix or second rigid-body integrator is introduced.
Physical impulses and velocities remain distinct from transformed determinacy
evidence; velocity-only admission retains the common-support-work restriction.

`resolve_articulated_elliptic_contacts_with` borrows original histories and
queries each geometrically active contact once, in canonical contact-ID order,
using actual incoming point velocities, normalized frames and both physical
witness positions. A plane source point is the normal projection, not its
arbitrary defining anchor. Inactive histories remain untouched. Mapped source
candidates return only with the globally admitted mechanical response, without
requiring a clone of every accepted history.

The surrounding event owner publishes all returned system responses and active
histories together. Source failure, cancellation, numerical refusal and late
physical audit failure return no partial candidate. This remains a fixed-pose
instantaneous impact law: no collision discovery, time advancement, pressure
from impulse or thermal allocation is inferred. The caller owns normal response
and source-frame transport. One-axis free slip uses the same native map and
complete original momentum/work checks; unproved determinacy still refuses.
Focused rigid roller-chain and off-centre articulated regressions check free-axis
momentum, shared normal/traction recoil, resisting-axis stick, permutation and
physical angular/work balances.

## Directional friction through sphere impact events

`sphere_plane::directional::impact` adds native and source-callback entry
points with an explicit instantaneous impulse ellipse. It reuses the original
normal restitution solve, contact detection, body impulse application,
prescribed-plane work and complete trajectory momentum/energy admission.
The existing sustained-only entry points retain their original impact law.

The native entry point takes separate sustained and impact laws/histories.
The impact law requires equal static/kinetic principal coefficients and
borrows its retained texture history without advancing slip or temperature.
The generic entry point exposes actual preimpact body state, event time,
normal velocities, admitted normal impulse and tangent mobility. The source
owner explicitly chooses instantaneous sampling; no duration, contact
pressure or finite travel is manufactured for a collision.

Support and impact callbacks share one private candidate history. A later
source refusal, bad impulse work or cancellation cannot publish preceding
private impacts or supported motion. Actual point velocity and midpoint
impulse loss are independently checked before the original body application
and whole-trajectory audit. Sustained heat and directly accumulated impact
loss remain separate outputs.

Four integration groups cover analytic impulse/spin/pose/loss and isotropic
recovery, detected impact followed by sustained arrest without invented
travel, rotated texture/moving-plane work, and late failure/cancellation/
forged-work refusal with deterministic retry. The single-sphere prescribed-
plane geometry admits decoupled normal and tangent mobility; this entry point
does not extend the contact law to general coupled rigid-body impacts.

## Finite Hertz shear/spin through sphere impact trajectories

`sphere_plane::hertz_impacts` applies the original finite-patch Hertz impulse
law at actual detected sphere/plane impacts. The real mass, isotropic inertia,
sphere radius and contact arm determine the two shear mobilities and normal
spin mobility. Normal restitution remains with the original normal owner;
normal and shear/spin responses decouple for this particular geometry. The
native Hertz surface couples shear and spin in one impulse budget, retaining
pure spinning, mixed sliding/spinning, complete arrest and zero-friction limits.
No independent saturated torque is added to an already saturated shear impulse.

One native momentum update applies the point impulse and its free normal-axis
angular impulse together. The original event receipt exposes that free couple,
and the whole-trajectory angular balance includes it about the world origin.
Changed spin propagates through the actual subsequent attitude, flight and
support motion. The independent unscaled midpoint work check includes both
shear and spin, before the original native endpoint and full-trajectory
geometry/momentum/energy gates. A translating plane supplies point-impulse
work; its zero angular velocity supplies no free-couple work. Held force and
motor work remain separate. Directly accumulated impact loss includes normal
plus shear/spin loss and is already part of the trajectory's total dissipation.
It is not an extra heat allocation.

The effective Hertz radius and normalized pressure-impulse profile are supplied
physical assumptions, never inferred from undeformed rigid sphere geometry.
The source-callback entry point receives actual incoming slip/spin, event time,
normal impulse and native body state, and returns an original admitted impulse
law plus candidate source history. It does not supply an arbitrary impulse.
There is no invented impact duration, normal pressure, finite travel or thermal
partition. Persistent support retains the separately supplied point-Coulomb
static/kinetic law and its explicitly zero torsional resistance; the impact
footprint is not extrapolated into a sustained constitutive model.

Only complete trajectory acceptance publishes source history and motion.
Later source failure, native law refusal, exhausted event/impulse budgets or
cancellation discard earlier private candidates. Six regression groups cover
analytic spin braking/arrest and attitude, mixed-budget recoil, detected
motor-driven impacts and continuation, rotation/boost/plane-work covariance,
separate sustained support, event-dependent source radius and atomic refusal.

## Finite Hertz shear/spin during persistent sphere support

`sphere_plane::hertz_support` advances the original sphere pose, linear and
angular momenta under one distributed circular Hertz static/kinetic wrench
surface. Shear and normal-axis spin share that surface; independent saturated
force/torque limits and a fitted wrench ellipse are not substituted. The
original support reaction supplies the load, and the radius is a separately
declared physical input. A source-returned `HertzWrenchContact` must retain that
actual load. The scene's point-impact model remains independently selected;
this support entry does not extrapolate its sustained patch into an impact.

The free spin couple enters the same original constant-wrench segment as its
tangential force. Actual spin drives attitude, and the free angular impulse
enters the whole world-origin angular balance. Actual contact-relative shear
and normal spin determine friction work. A translating plane supplies its real
force/velocity work and has no invented torsional work. External motor work,
physical sustained loss and original impact loss retain separate accounts.

Pure spin and collinear pure shear use the original exact Hertz capacities and
analytic first-zero times. Static holding and breakaway test the complete force
and moment demand against the original joint static surface. The accepted
`SpherePlaneHertzState` distinguishes resolved zero motion from durable static
capture: an overdriven zero is continued into reversal, including an interval
ending exactly at that zero, without resetting native rounded momenta. Generic
mixed motion uses bounded velocity-implicit projection on the original kinetic
surface, with first-order continuum convergence and timestep-boundary capture.
Moving capture additionally requires interior kinetic-demand feasibility; a
coarse step cannot capture within the static/kinetic hysteresis band. Trials
whose held wrench assists the incoming motion refine or refuse. The actual
linear-velocity segment's midpoint work owns dissipation; the constitutive
projection's endpoint work is not a second heat source.

Retained roundoff at static capture can have small signed native constraint
work. `SpherePlaneAdvance::static_contact_work_residual` reports this account
explicitly, and the original global energy gate closes with it separately from
nonnegative physical loss. Each static segment and the sum of absolute static
defects have independent caller bounds. No momentum projection, clipped heat
or relaxed global conservation tolerance hides this defect. Original point and
material rolling paths retain exactly zero in the new account.

Source phases, mixed-motion steps and attitude steps have distinct controls.
The source callback receives the actual body, point, frame, load, shear/spin,
driving acceleration and interval; it returns an original law plus private
candidate history. Only complete original geometry, momentum, work and budget
acceptance publishes motion and history. Refining the mechanical timestep does
not resolve a changing held source. Hot-step work uses fixed arrays; support
transitions and source work are bounded. This model adds no retained partial
slip, deformation rolling loss, normal compliance, temperature evolution,
contact-radius evolution or heat partition.

Eight native regression groups cover analytic spin stopping/attitude, static
hysteresis and cumulative static-work refusal; opposed-motor zero/reversal and
split continuation; recovery of the original pure-shear trajectory; joint
static holding/breakaway; first-order refinement against independent continuous
RK4 COM/spin/position/dissipation flow; actual source-dependent radius and atomic
refusals; detected impact-to-support and lift-off; and arbitrary rotation,
shared drift and real moving-plane work. Together with the existing native
sphere-plane suites, 71 tests pass. These are bounded model checks, not physical
calibration, full-workspace Cargo or performance qualification.

## One finite-patch trajectory through impact and sustained support

`sphere_plane::hertz_contact` composes the original Hertz impulse event and
sustained Hertz shear/spin driver inside the same native event loop. Impact
and support retain separately declared radii, coefficients and dimensioned
controls. Neither law is inferred from the other. The original normal impulse,
restitution and event timing remain unchanged; the actual angular and linear
momenta after each impact start the subsequent support phase. No second
impulse, body reconstruction or external event-splitting clock joins them.

The original impact phase, solve and independent unscaled endpoint/work audit
are shared with `hertz_impacts`. The original support candidate, capture/reversal
logic, static-work account and whole-motion geometry/momentum/energy gates are
shared with `hertz_support`. Mixed support still requires first-order time
refinement. Impact receipts retain their original pressure-impulse assumption;
support retains its independently supplied load/radius and force law.

`advance_sphere_plane_hertz_contact_with` stages two independent source-history
types and preserves their original callback/error types. Actual event ordering
determines which owner is queried. Both histories are cloned once, after control
admission. The complete trajectory and final cancellation check must succeed
before either caller history changes; both swaps finish before a displaced
history can run its destructor. A later support, impact, source or native
refusal therefore cannot publish an earlier event's history. Caller-owned
callback effects retain their original boundary.

The returned `support.motion` is the single accepted native trajectory. Its
physical loss already includes both `support.impact_loss_j` and
`support.sustained_friction_loss_j`; event receipts are not added again as heat.
The composition adds no heat partition, finite impact duration, normal
compliance, measured calibration or new constitutive approximation.

Three G1/G3/G4 integration regressions verify analytic flight, detected impact,
support spin arrest and attitude with distinct radii and 6.6 J total physical
loss; inclusive-event split continuation without a duplicate impact; actual
post-impact source kinematics and both-history rollback after late support
refusal/cancellation; and a second-bounce source refusal followed by exact
retry. All 74 tests in the complete focused native sphere-plane target pass,
including the original impact and sustained-support regressions.

## Directional friction on moving sphere/plane trajectories

`sphere_plane::directional` advances the original sphere's position,
attitude and world momenta through the native directional law's accepted
held-traction segments. Its native and source-callback entry points reuse
the original sphere/plane engine for normal support, flight, impact,
lift-off/recontact, supported-path geometry and whole-trajectory momentum
and energy checks. Normal reaction follows the applied force; changing that
force between calls can release or restore contact.

The plane owns a fixed texture frame. Each supported phase supplies its
actual contact velocity, sphere mobility, driving acceleration and normal
reaction to the constitutive owner. The callback resolves source trials
privately and replays only accepted native segments. Source-phase duration,
native curved-sliding resolution and attitude resolution have independent
controls; refining one does not establish accuracy of the others. Constant-
direction intervals retain native analytic arrest, while curved sliding
requires timestep refinement.

Actual point kinematics, integrated impulses, displacement and traction-slip
work are checked against every accepted segment and the aggregate report.
Only complete trajectory acceptance publishes the returned body state and
caller-owned constitutive history. Source failure, exhausted budgets or
cancellation leave the caller's history unchanged. Sustained friction heat
and instantaneous impact loss have separate directly accumulated outputs;
external-force and prescribed-plane work remain with the original engine.

The ideal sphere has no inferred contact area or pressure. A source callback
may obtain temperature and pressure from declared external owners; this
module does not assign heat to thermal cells. Plane texture does not rotate
with sphere attitude. Impact friction retains the original scene law, and
normal spin has no torsional resistance. This is a single sphere against a
prescribed plane, with constant forces and plane velocity within each call.

Six moving-sphere groups cover analytic arrest/pose/work, static hold and
breakaway, covariance, curved-path refinement, flight/impact loss separation,
and atomic failure/retry. Cross-crate source tests exercise speed-dependent
material claims through actual motion, lift-off/recontact and late source
failure without changing the source/history owner.

## Explicit drive in directional rigid contact

`tangential_friction::directional::driven` applies held external world forces
and body-frame torques together with the original integrated directional
traction to the original finite body pair. It reuses full off-centre tensor
mobility and point kinematics. The native and source-callback entry points
publish body momenta, accepted friction history and Resting/Sliding mode
together after independent impulse, endpoint and energy checks.

External work is evaluated through the mass-metric decomposition: free
full-body work minus free relative-contact work plus the driven trajectory's
actual relative work. This retains common-body and complementary rotational
motion, so static relative hold can gain kinetic energy without producing
frictional heat. The independent balance is `delta K + D - W_external`;
native effective drive, body momentum changes and callback work are audited
separately. A force pre-kick is not substituted for the simultaneous drive.

Four driven rigid groups cover common-body static work, breakaway, full-tensor
off-centre response, frame covariance, zero-drive recovery, and rejection of
false work/drive/rest receipts or cancellation without publication. Poses,
inertias, contact arms and frame are held during this phase. Actual normal
recoil is reported; normal constraints, pose evolution and material/thermal
refresh retain their existing owners.

## Source callbacks for continuous directional trajectories

`apply_rigid_directional_trajectory_with` supplies the original pair's full
tensor mobility and point velocity to a constitutive source callback. The
callback returns the native integrated trajectory plus its source/history
candidate. The native constant-law entry point delegates to this same owner;
there is one body mobility calculation, impulse application and audit path.

Only the original two bodies are staged. Opposing integrated impulses, actual
endpoint velocity, individual and total momenta, impulse/kinetic work and
nonnegative physical heat are checked before those bodies and the callback
history publish together. A source refusal, bad candidate or cancellation leaves
all accepted owners unchanged. The callback owns constitutive/source admission;
this adapter does not manufacture coefficients, resolve material names or
replace native curve work with average force times net travel.

Poses, arms, inertia and normal load are frozen. Normal recoil is reported, not
silently constrained. Changing sources and thermal storage remain external
owners. G1/G3/G4 native tests retain analytic and full-tensor arrest/continuation
checks and add a source-history callback, audit rejection and atomic retry.

## Purpose and layer

Layer: **L3 FLUX**. `fs-mbd` provides safe-Rust rigid-body dynamics in two
deliberately explicit lanes. The original lane covers one smooth,
unconstrained rigid body whose reference point is its centre of mass. The
`articulated` lane adds validated general spatial inertia, scalar joints, a
root-first rigid-body tree, Lie-group forward kinematics, recursive
Newton-Euler inverse dynamics, and Featherstone articulated-body forward
dynamics in linear time and storage.

The `robot_models` catalog is a provenance-bound construction layer over that
articulated owner. It currently transcribes two pinned upstream descriptions:
Unitree's current mode-11 G1 description with all 29 actuated joints, and the
seven-axis KUKA LBR iiwa 7 R800 description from `iiwa_stack`. It adds no
second pose, twist, wrench, inertia, joint, or articulated-model type.

Canonical poses, twists, wrenches, adjoints, and coadjoints remain owned by
`fs-ga`; the articulated lane consumes those types rather than creating a
parallel robot-math representation. The articulated owner now includes both a
prescribed-base boundary and an unconstrained free-flight boundary that solves
six unactuated base accelerations together with the scalar joint accelerations.
The dedicated `articulated::motion` owner advances continuous unconstrained
free-flight trajectories from those native dynamics. Contact, impacts and loop
constraints remain separate owner boundaries. Control policy is likewise a
separate owner: this crate consumes declared efforts and external wrenches but
invents neither.

The `goddard` module also exposes a narrow museum-facing composition for US
1,102,653. It reuses the single-body owner for torque-free primary and
gyroscope display poses, then reports only source sequence predicates and the
Claim 2 tapered-tube ratio. The older `step_goddard_rocket` liquid-nozzle
calculation remains an explicitly adjacent ideal calorically-perfect model;
it is not a model of the mechanism claimed in US 1,102,653. Its chamber
stagnation pressure and temperature, specific gas constant, constant cp/cv,
throat area, exit/throat area ratio and ambient pressure are all explicit.
The sonic throat determines mass flow; there is no independently prescribed
inconsistent fuel-flow input or temperature inferred from pressure.
An allocation-free, at-most-80-iteration inversion selects the supersonic
area-Mach branch (exactly sonic for unit area ratio). The exit obeys the same
ideal-gas energy and mass balances as the throat. Total thrust includes both
`mdot*ve` and `(pe-pa)*Ae`; specific impulse uses total thrust and standard
gravity, and no physical output is rounded for presentation. Deterministic
elementary functions and cancellation-safe near-sonic and gamma-near-one
expressions are used. Accepted gamma is `(1,5/3]`, area ratio is `[1,1000]`,
and ambient pressure must lie in `[0,pe]`; overexpanded operation refuses
until a shock/separation model is available. Invalid or unrepresentable states
return `GoddardRocketError`, never a partial result. The G1 tests use the
analytic gamma=7/5, Mach=2 solution, the sonic and gamma-to-one limits,
mass/enthalpy balance, pressure-thrust differences and geometric scaling.

The `daimler` module exposes a second narrow museum-facing composition for US
361,931. It constructs the printed longitudinal propeller-shaft freedom through
the generic articulated prismatic-joint owner, then reports the mutually
exclusive ahead, neutral, and astern contact topology plus the passive/pumped
cooling-path selection. Because the grant prints no travel, load, friction,
speed, flow, or power values, this composition publishes normalized topology
only and no quantitative performance.

The `davinci` module exposes the source-bounded joint topology printed in Figs.
2 and 2A of US 6,331,181. It composes five generic revolute joints and one
normalized prismatic insertion joint, retains the claimed compatibility
identifier as a predicate, and deliberately publishes no link dimensions,
motor data, force, power, accuracy, or clinical-performance quantity.

The `howe` module exposes the source-order joint topology printed in US 4,750:
one main shaft constrains the curved needle arm, picker-driven shuttle in its
horizontal trough, loop-lifting rod, and baster-plate feed to one prescribed
drive coordinate. It publishes normalized kinematics, interlock predicates,
and only the two local dimensions printed by the grant; it publishes no
invented machine dimension, force, torque, speed, friction, or power.

The `planar_drive` module owns a reusable constant-twist SE(2) update for
prescribed left/right wheel speeds. The `roomba` module composes that generic
owner with the intersecting emitter/detector field and change-of-direction
logic of US 6,594,844. Its room and low-furniture boundaries perform bounded
kinematic non-penetration projection for a museum display; they do not add a
contact-force, impact, tire, friction, traction, or cleaning-performance model.

The `salisbury` module exposes the source-bounded joint topology and static
tendon law printed in US 4,921,293. It composes three generic revolute joints
for each of three digits and evaluates the grant's three Figure 3 torque
equations from caller-declared SI tensions and an explicitly illustrative
pulley-radius scale. It publishes neither a fabricated historical dimension nor
a dynamic, contact, grasp-stability, or force-closure result.

## Public types and semantics

- `Vec3` is a small Cartesian vector. Every public field that carries a vector
  names its frame in the field name or documentation.
- `UnitQuaternion` is a normalized body-to-world quaternion. Construction uses
  scaled normalization, so finite huge and subnormal components do not overflow
  or underflow the normalization intermediate; it rejects only zero/non-finite
  input. The representative has the first nonzero component of `(w, x, y, z)`
  positive, making the double cover deterministic.
- `Pose::new` combines a finite world-space centre-of-mass position with that
  orientation. Its fields are private; `position_world()` and `orientation()`
  expose copies of the validated values.
- `MassProperties::new(mass, center_of_mass_body, principal_inertia_body)`
  accepts finite positive mass and principal moments satisfying the rigid-body
  triangle inequalities. The centre of mass must be exactly zero in the
  principal frame; offset-reference spatial inertia is refused.
- `MassProperties::from_principal_frame` admits known principal moments and
  their fixed orientation in the actual body coordinates.
  `from_inertia_tensor` uses the shared `fs_time::lie::gyrostat::InertiaTensor3` admission
  once. Both retain the same rigid-body triangle and COM-reference rules.
  `principal_inertia_body()` returns eigenvalues; `principal_axes_to_body()`
  supplies their axes and `inertia_tensor_body()` supplies the full tensor.
- `RigidBodyState` stores world-frame linear momentum and body-frame angular
  momentum. Its fields are private and exposed through frame-named getters; it
  has no independent velocity or angular-velocity state.
- `Pose::point_world_from_body` and `Pose::point_body_from_world` are checked
  centre-of-mass-relative point transforms. `UnitQuaternion` provides the
  corresponding canonical world-to-body vector transform.
- `RigidBodyState::point_kinematics` returns a validated `PointKinematics`
  record with the body/world contact arm, point position, centre-of-mass
  velocity, angular velocity, and material-point velocity
  `v_com + omega_world cross r_world`.
- `RigidBodyState::directional_effective_mass` returns a normalized-direction
  `DirectionalEffectiveMass`. It is the unconstrained scalar free-body response
  at the declared body point; zero/non-finite directions and non-positive or
  non-finite derived denominators refuse.
- `RigidBodyState::apply_impulse_at_body_point` applies one world-frame impulse
  at a centre-of-mass-relative body arm. It updates `p_world += J_world` and
  `L_body += r_body cross J_body` without moving the pose, and returns a full
  `ImpulseReceipt` with before/after point kinematics and kinetic/work
  diagnostics. `apply_force_at_body_point` explicitly converts a finite force
  held for a positive declared duration into that same atomic momentum event;
  it intentionally does not pretend to integrate pose evolution during the
  duration.
- `apply_equal_and_opposite_impulse_at_body_points` returns two atomic event
  receipts plus algebraic impulse balance, floating-point linear-momentum
  change, and paired kinetic/work accounting. It does not infer that the two
  declared points coincide or constitute a physical contact.
- `Wrench` carries a constant world-frame force and declared-body-frame torque
  for one step. `Gravity` is a validated uniform world-frame acceleration with
  a private field and `acceleration_world()` getter.
- `RigidBodyIntegrator::step` returns a `StepReceipt` with the before/after
  states and diagnostics. Translation uses the midpoint linear momentum under
  the constant total force. Body angular momentum solves implicit midpoint for
  `Ldot = L × I⁻¹L + tau`, with the actual endpoint mean as its conjugate
  momentum. A matching right-composed Cayley attitude update preserves world
  angular momentum in torque-free asymmetric rotation. The quadratic kinetic
  energy balance is `K1 - K0 = dt * tau dot I^-1((L0 + L1)/2)`, up to the
  nonlinear-solve tolerance and floating-point roundoff. Spherical bodies and
  single-principal-axis rotations use the direct momentum update and quaternion
  exponential; constant-axis spin therefore remains exact at finite elapsed
  time rather than acquiring a Cayley phase error.
- `RigidBodyIntegrator::diagnostics` returns
  `Result<DynamicsDiagnostics, DynamicsError>`. It revalidates its state,
  gravity, and mass properties and refuses non-finite derived energies or world
  angular momentum rather than publishing an `inf`/`NaN` ledger.
- `RigidBodyIntegrator::step_symplectic` is an explicit second-order alternative
  to the energy/momentum midpoint lane. Body-torque half-kicks surround a
  symmetric Hamiltonian free-rotation split, preserving spatial angular momentum
  to roundoff without a nonlinear solve. A commuting spherical term uses the
  median principal inertia, giving exact spherical and axisymmetric free motion;
  asymmetric bodies retain a bounded second-order energy defect at adequately
  resolved fixed steps. It does not inherit `step`'s exact energy/work identity.
- `RigidBodyIntegrator::step_kick_drift` applies the full declared torque kick,
  then uses the same free-rotation map as `step_symplectic`. Its forced
  composition remains first order; torque-free rotation is second order and
  preserves spatial momentum. Each body torque kick contributes its world
  impulse at that kick's actual attitude.
- `RigidBodyIntegrator::advance` invokes its cancellation callback before each
  whole step. It returns `AdvanceOutcome::Cancelled` with the last fully
  committed state and diagnostics when cancellation is observed.
- `goddard::step_goddard_apparatus` admits declared primary/gyroscope RPM,
  elapsed seconds, the exact `L/D` source ratio, auxiliary-release progress,
  the source's "substantially consumed" state, and gyroscope presence. It
  returns normalized torque-free `fs-mbd` quaternions, angular velocities,
  Claim 1/2 sequence predicates, and ideal instrument-support isolation only
  while the declared gyroscope is both present and spinning. A present but
  stopped gyroscope shares the primary body's world rate. The normalized
  spherical mass properties affect no published force, energy, or trajectory
  quantity.
- `daimler::step_daimler_marine` admits exactly three shaft selections
  (`-1` astern, `0` neutral, `1` ahead) and optional pump state. It constructs
  a normalized one-DoF `articulated::JointModel::prismatic` along the vessel's
  longitudinal axis. Ahead returns negative sternward translation (movement
  toward the motor), closes only `a/a²`, and retains the source-stated thrust
  contact predicate; astern opens that contact and closes only the
  `e¹/e²`-with-`a²/c` path. The motor sign remains positive in every state,
  matching the continuously one-direction source statement. Pump activation
  never removes the printed fore/aft pipe path.
- `davinci::step_davinci_topology` admits finite support/tool angles, one
  dimensionless insertion coordinate in `[-1, 1]`, and the tool compatibility
  identifier predicate. Its revolute axes and prismatic tool axis come from
  generic `articulated::JointModel` owners. Compatibility is reported rather
  than reinterpreted as a fabricated motor interlock, and normalized insertion
  is not represented as measured travel.
- `howe::step_howe_topology` admits a finite main-shaft angle, normalized
  displayed loop slack in `[0, 1]`, and the Claim 1 combination predicate. It
  composes seven generic scalar joint coordinates constrained to one prescribed
  drive coordinate, returns the source-order needle/loop/shuttle/feed phase,
  and refuses to label declared display slack as a printed dimension.
- `planar_drive::step_differential_drive` admits a finite planar axle-midpoint
  pose, finite prescribed left/right tangential wheel speeds, positive track
  width and wheel radius, and a fixed step in `(0, 0.25]` seconds. It applies
  the exact constant-twist SE(2) exponential and updates wheel display angles;
  equal speeds translate, opposite speeds spin about the axle midpoint.
- `roomba::step_roomba` admits the shared room dimensions, procedural low-solid
  footprints, reader-selected speed/turn controls, optical sensor geometry,
  and the previous deterministic tape state. It composes the generic planar
  drive, reports the source-bounded surface-field overlap and wall-field
  redirect, projects a circular display bumper outside room/low-solid
  boundaries, and returns stable contact indices (`-1` clear, `-2` room wall,
  otherwise the caller's collider index). At most 64 colliders are admitted.
- `salisbury::step_salisbury_hand` admits four finite non-negative cable-end
  tensions in newtons, a finite positive visitor-declared R2 scale in metres,
  and a Claim 2 first-idler predicate. It composes nine scalar revolute joint
  coordinates through `articulated::JointModel`, reports twelve cable ends,
  and emits the root-first parent map
  `[-1, 0, 1, -1, 3, 4, -1, 6, 7]`: each `-1` is one digit's palm anchor,
  and each following pair is its serial Axis-2/Axis-3 chain. This topology
  proves attachment without fabricating a link length, mass, or inertia. It
  evaluates the exact printed relations
  `Torque1 = -T1 R1 + T2 R2 + T3 R2 - T4 R1`,
  `Torque2 = T1 R3 + T2 R2 - T3 R2 - T4 R3`, and
  `Torque3 = T2 R2 - T3 R2`. The study ratios `R1=1.2 R2` and `R3=1.4 R2`
  only retain Figure 3's depicted ordering; they are not historical dimensions.
- `articulated::SpatialInertia` accepts positive mass, a finite centre-of-mass
  offset, and a full symmetric centre-of-mass inertia. It validates positive
  definiteness and the principal-moment triangle inequalities before exposing
  the corresponding 6x6 spatial inertia, momentum, or kinetic energy.
- `articulated::JointModel` has validated revolute, prismatic, and helical
  constructors plus a fixed-joint constant. Its axis and motion subspace cannot
  be forged through public fields. `JointLimits` carries inclusive position and
  symmetric speed/actuator-effort bounds.
- `articulated::ArticulatedModel` requires one root at link zero, root-first
  topological ordering, preceding parents, unique nonempty names, and compact
  scalar-DoF indexing. The prescribed-base API describes the frame above the
  root with `BaseState`, including its body twist and non-gravity body
  acceleration; that boundary does not solve the base motion.
- `articulated::FreeFloatingBaseState` carries the canonical `fs-ga` `Se3`
  world pose and body-coordinate `Twist` of the frame above the root. It has no
  acceleration field because `free_floating_forward_dynamics` solves the six
  unactuated base coordinates. A free-floating model must have a fixed root
  joint; a scalar root joint would redundantly parameterize motion already
  represented by the base pose and is refused.
- `articulated::forward_kinematics` returns world-from-link `Se3` poses and
  body-coordinate twists. `inverse_dynamics` implements a recursive
  Newton-Euler pass and reports required generalized effort. `forward_dynamics`
  implements Featherstone's articulated-body algorithm without constructing a
  dense generalized mass matrix; supplied actuator efforts are checked against
  declared limits.
- `articulated::free_floating_forward_dynamics` reuses the same linear-time
  Featherstone factor/back-substitution passes, accumulates one 6×6 articulated
  inertia at the root, and solves that fixed-size system for the physical base
  spatial acceleration. Uniform world gravity is applied as physical link body
  wrenches, so `FreeFloatingForwardDynamics` reports Featherstone spatial
  accelerations with physical gravity included rather than the prescribed API's
  gravity-shifted recursion convention. These are explicitly named
  `base_spatial_acceleration_body` and `body_spatial_acceleration`: for body
  twist `[omega, v]`, the ordinary Cartesian acceleration of the frame origin
  is `a.linear + omega x v`, exposed by
  `origin_linear_acceleration_body`. External wrenches remain body-coordinate
  per-link values and the returned scalar accelerations retain compact joint
  order.
- `ArticulatedModel::free_floating_complexity` admits only the same fixed-root
  topology as the solver, then records six base DoFs, one fixed 36-entry root
  solve, zero dense generalized-matrix entries, and the shared linear tree
  working set.
- `robot_models::unitree_g1_29dof` builds a fixed-pelvis tree in the current
  mode-11 source order: six left-leg, six right-leg, three waist, seven
  left-arm, and seven right-arm DoFs. It uses all 29 actuated link inertias,
  joint origins, axes, and hard limits from the official Unitree
  `g1_29dof_mode_11.urdf` pinned in the returned metadata. The pinned Unitree
  README identifies mode 11 as the up-to-date 29-DoF configuration.
- `robot_models::kuka_lbr_iiwa7_r800` builds the fixed-base seven-axis chain
  from `iiwa_stack`'s pinned `iiwa7.xacro`. Its 300 N m effort and 10 rad/s
  velocity limits are the Xacro's generic macro defaults, not claimed KUKA
  hardware limits.
- `robot_models::CatalogRobotModel` binds a validated `ArticulatedModel` to its
  stable compact joint order and immutable `RobotModelMetadata`. The metadata
  records exact pinned URLs, revisions, Git blob identities, source status,
  units, derivation, and material omissions. `ROBOT_MODEL_CATALOG_VERSION`
  versions this typed in-source layout; it is not a serialized URDF schema.
- `robot_models::G1ResidualPolicy` admits the catalog-owned 5,040-weight
  actuator-major policy once, then evaluates the exact 42-signal by eight-phase
  feature map for each observation without repeatedly validating immutable
  weights. Deterministic `tanh` bounds every normalized actuator residual.
- `robot_models::g1_policy_phase_basis` publicly owns the exact periodic basis
  `[1, sin(phi), cos(phi), sin(2phi), cos(2phi), sin(3phi), cos(3phi),
  sin(4phi)]`. Features are signal-major within each 336-weight actuator row;
  curriculum code selects owner coordinates from this layout rather than
  reproducing a second robot-policy representation.

## Work-conjugate finite steps for held spatial forces

`RigidBodyIntegrator::step_spatial_midpoint` consumes the same
`spatial_loads::PreparedSpatialLoads` as the existing Verlet lane. Held material
point forces now have an explicitly selected finite-step lane whose native
torque work equals force dotted with the actual material-point displacement.
For `Rbar = (R0+R1)/2`, the mean spatial moment uses the mean material arms and
is pulled back by `Rbar^-1`. A normalized midpoint attitude would change that
work and its spatial angular impulse. World-held and body-held free couples
retain their separate frames; their work uses the conjugate Cayley increment
and remains a second-order continuous-work quadrature.

A stack-only three-variable Newton solve uses an analytic full-inertia/load
Jacobian and the original three prepared first moments. Per-step work and
storage are independent of the point-force count. The selected torque is
replayed through the original `step_cayley` owner; its actual returned momentum
must satisfy the same held-load equation before publication. No momentum,
energy or pose projection is added. Differential quadratic-energy evaluation
avoids subtracting two large total kinetic energies.

Explicit limits bound Newton corrections, backtracking, angular-momentum
residual and the Cayley increment. Cancellation at every solver trial and
before final return, failed native replay, nonfinite arithmetic or exhausted
limits return no updated state. The full-inertia work balance closes to solve
tolerance and roundoff; phase/trajectory accuracy still requires refinement.
This is neither a symplecticity/global-uniqueness certificate nor a contact or
state-dependent-force model, and numerical defects are never assigned as heat.

Five G1/G3/G4 groups in `tests/spatial_loads/midpoint.rs` cover analytic
translation/axial torque, physical point-force pendulum work and second-order
refinement, long-time energy, full-tensor reversibility and coordinate
covariance, spatial impulses, original held-body-load equality, cancellation at
every publication boundary and bounded/nonfinite refusal with deterministic retry.

## Generalized articulated coupling

`articulated::generalized::free_floating_terms` exposes the original admitted
tree's complete free-base mass matrix and velocity bias in
`[omega_base, v_base, joint_rate]` coordinates. Its equation is
`M(q) udot + bias(q,u) = applied_force`. The base linear derivative is in body
coordinates; Cartesian origin acceleration adds `omega cross v` before world
rotation. Gravity is excluded from bias and explicitly projected at each link
COM by `project_uniform_gravity`.

The native composite rigid-body pass retains all locked-base, base/joint and
joint/joint mass terms, including moving COM offsets and rotated full link
inertias. The original Newton-Euler recursion supplies bias and each link's
zero-generalized-acceleration term. `GeneralizedKinematics` shares actual link
poses, twists, analytic Jacobian columns and these acceleration terms with
additional fluid/structural kinetic operators. It also projects body-held or
world-axis link-origin wrenches by virtual work. No finite differences, unit-rate
admission probes, duplicated Newton-Euler engine or additional time integrator
are introduced.

This explicit coupling seam admits at most 4096 links and 64 scalar joints;
native unconstrained ABA retains its existing linear-storage path. Projection
work/storage are bounded by `O(links*(6+joints) + (6+joints)^2)` and cancellation
is polled around the original indivisible FK/RNEA passes and during projection.
Native fixed-root, physical-inertia, configuration, rate and finite-output
admissions remain in force. `articulated_generalized` verifies independent
serial-rotor and moving-slider mass/recoil formulas, mixed-joint energy,
CRBA/RNEA/ABA agreement, distributed gravity, load-frame covariance, caps and
atomic cancellation. Extraction itself makes no time-integration conservation
claim; its consumer owns that numerical approximation.

## Rigid dynamics in arbitrary body coordinates

### Spatial torques and forces on moving material points

`spatial_loads::PreparedSpatialLoads` combines the original `Wrench` (world
COM force and body-held free couple), an independently world-held free couple,
and world forces at fixed body material points. The latter produce
`sum (R r) cross F`; their moment is reevaluated at the actual attitude.
The declared body coordinates, original full inertia and center-of-mass
reference are shared with `RigidBodyState`, not replaced by a new state owner.

Preparation is a deterministic compensated O(points) reduction to the force
resultant and three first-moment vectors, with constant storage and no heap
allocation. At most 16,384 points are admitted; cancellation is checked before
work, every 64 points and before return. All input vectors, moment products and
reductions must be finite and representable. A compiled load can be reused at
any body attitude; both moment evaluation and each subsequent step have
constant work independent of the number of input points. This is an algorithmic
complexity statement, not a machine-specific performance measurement.

First moments retain both multiplication and summation tails using the existing
double-double arithmetic owner. Torque contractions, midpoint load/Jacobian
evaluation and endpoint force work carry those tails until forming the native
resultant. Thus nearly radial forces and large opposing point loads retain a
representable residual couple instead of silently losing its angular impulse.
Exactly cancelling colocated forces still produce zero force and torque. This
does not add precision to the stored native state or certify arbitrary dynamic
range; nonfinite products and reductions retain their original refusals.

`RigidBodyIntegrator::step_spatial` delegates its endpoint momentum half-kicks
and Hamiltonian free-rotation map to the original `step_verlet_with` owner.
It applies gravity once at the COM
and integrates the constant resultant force exactly. A pure world-held couple
obeys `R1 L1 - R0 L0 = dt tau_world` to roundoff; body couples and point forces
contribute their endpoint spatial moments. The receipt reports actual linear
and total angular momentum change minus that applied impulse, including the
orbital moment of force about the world origin.

Held point-force work is reconstructed from actual endpoint point displacement,
using the same first moments. Free-couple work uses the two kick power
quadratures, with the original full inverse inertia. Mechanical-energy change
minus these work terms remains an explicit numerical defect, never heat or
physical damping. The held-load composition is second order; exact forced
energy preservation is not claimed. Conservative point forces retain the
Hamiltonian split. State-dependent drag, springs, contacts, changing fields or
deforming attachments require their own stage evaluation and are not modeled
by freezing a prepared load. `advance_spatial` commits only complete accepted
steps and checks cancellation at each whole-step boundary.

Focused G0/G1/G3 regressions in `tests/spatial_loads.rs` exercise world angular
impulse, exact axial forcing, direct point-moment/work reconstruction,
full-inertia coordinate changes, refinement and atomic refusal/cancellation.

### Shared full-inertia representation

The native `MassProperties` and `RigidBodyState` owners now accept fixed full
COM inertia tensors while keeping geometry, contact arms, body torques and
angular momentum in the caller's declared body coordinates. Known principal
frames require no eigensolve; a full matrix is decomposed once by the existing
bounded shared owner, with its numerical reconstruction gate. Positive kinetic
definiteness alone is insufficient here: the physical principal-moment triangle
inequalities also apply. The cached representation adds no per-step allocation,
factorization, or eigensolve. The old diagonal constructor and its integration
paths retain their original numerical results.

Both implicit midpoint/Cayley and symplectic free-rotation maps reuse their
original principal-frame numerical owners through a fixed change of basis.
Their attitude increments, momentum and torque use the same transformation.
Spherical/axisymmetric shortcuts therefore refer to actual principal directions,
not arbitrary geometry axes. Point kinematics, kinetic energy, directional
effective mass and impulse work all use the same full inverse tensor.
Normal/Coulomb, finite-patch shear/torsion, Hertz sliding/spin and constrained
rolling mobility
assemble rotational Gram factors `diag(1/sqrt(I_principal)) Q^T`; products of
inertia cannot disappear merely because a geometry frame was relabelled.

`angular_momentum_body_checked` supplies the matching velocity-to-momentum
conversion for callers; componentwise multiplication by the returned principal
eigenvalues is valid only in the associated principal frame. Focused native
tests compare the same physical trajectories and contacts under different body
coordinate choices, independently reconstruct off-centre impulse response,
check the nonprincipal symmetric-top closed form, preserve diagonal results
bit for bit, and reject invalid/nonphysical tensors. Existing energy,
friction and rolling regressions remain applicable. This is fixed-inertia
rigid motion, not a changing-mass, moving-COM, deforming-body, interval-inertia,
material-calibration or arbitrary offset-reference spatial-inertia model.

## Invariants

- Public construction rejects non-finite mass-property, state, gravity, wrench,
  duration, and quaternion inputs at the relevant boundary. State, pose, and
  gravity fields cannot be forged externally; public integrator operations also
  revalidate their inputs before producing a receipt or diagnostic.
- The public orientation remains normalized and uses a canonical quaternion
  sign after every successful attitude update.
- The declared force and gravity are world-frame quantities; the declared
  torque and angular momentum are declared-body-frame quantities. Diagnostics
  rotate the latter to world space before adding orbital angular momentum.
- Each successful `step` constructs and returns a complete new state. In
  `advance`, only complete preceding steps are committed when cancellation is
  observed.
- Point kinematics, effective-mass queries, and event receipts validate all
  state, mass-property, arm, direction, force, impulse, and derived finite
  values before returning. The event APIs consume state by value and return a
  new state only inside a complete receipt, so a refusal leaves externally held
  state unchanged.
- Articulated model construction validates topology, names, axes, joint limits,
  and physical spatial inertia. Evaluation validates dimensions, finite values,
  joint positions and speeds, external wrenches, and forward-dynamics actuator
  efforts before running either recursion. Spatial-inertia matrices, momenta,
  kinetic energies, transformed twists/wrenches, articulated inertias, and
  returned accelerations are checked after arithmetic; large finite input that
  overflows a derived quantity refuses instead of publishing `NaN` or infinity.
- Free-flight evaluation additionally validates the canonical base twist,
  requires a fixed root joint, and transactionally refuses a non-symmetric,
  singular, non-finite, or ill-conditioned accumulated root inertia. It never
  mutates caller-owned state or publishes partial accelerations on refusal.
- Catalog builders deterministically preserve their declared root-first link
  topology and compact source-joint order. Their retained numeric records use
  URDF SI units and pass the same articulated topology, limit, pose, and
  physical-inertia validation as caller-built models.
- Every one of the 5,040 flat G1 policy coordinates is covered by an exhaustive
  owner test that perturbs it and verifies that exactly the corresponding
  actuator response changes under an observation with nonzero basis support.

## Error model

`DynamicsError` is the refusal channel for non-finite input or derived result,
unrepresentable finite magnitude, invalid mass, invalid or physically
inconsistent principal inertia, invalid orientation, invalid directional
effective mass, unsupported non-centre-of-mass reference point, and invalid
duration. `RotationalSolveFailed` reports a singular or unconverged implicit
rotational candidate and requests a smaller timestep. The fixed-size solve
normalizes momentum, uses an analytic 3x3 Jacobian with deterministic partial
pivoting, and permits at most 16 Newton iterations with at most eight trial
fractions each. It allocates no heap storage and never publishes an unconverged
state. There are no panics in production paths for those invalid inputs.
Arithmetic overflow or a non-finite attitude update, kinematic transform,
effective-mass query, or impulse-work result is returned through the checked
boundary; already completed earlier steps of `advance` remain committed.

The free-floating root solve scale-normalizes and symmetrizes the accumulated
6×6 physical inertia within a stated rounding tolerance, performs a private
fixed-size Cholesky solve, and computes a deterministic infinity-norm condition
estimate from six checked triangular solves. Non-positive/small pivots and
condition estimates above `1e12` are structured refusals. This fixed-size
helper is private and is not a second public linear-algebra API.

## Determinism class

The step order, frame transformations, quaternion sign convention, and
cancellation boundary are deterministic functions of the admitted inputs. The
implementation currently uses the Rust standard library's floating-point
`sqrt`, `sin`, and `cos`; therefore this contract makes no cross-ISA bitwise
claim and retains no cross-ISA evidence. Same-platform repeatability is an
implementation property to be exercised by consumers, not a retained G5
promotion from this eleven-test foundation.

## Cancellation behavior

This crate is synchronous and creates no threads, tasks, or I/O. `advance`
checks a caller-supplied cancellation predicate immediately before every
candidate step. A cancelled iteration does not start that step, so it cannot
partially change position, momentum, orientation, or its diagnostic result.
There is no `Cx` integration, deadline budget, allocation budget, drain phase,
or asynchronous cancellation protocol in this initial core.

## Unsafe boundary

None. The crate denies unsafe code and contains no unsafe blocks.

## Feature flags and dependencies

No feature flags. The articulated lane depends on `fs-ga` for its canonical
Lie-group and spatial-vector types; `robot_models` reuses that same dependency
and the articulated types. The legacy single-body lane otherwise continues to
depend only on `std`/`core`.

`fs-mbd` deliberately does not depend on `fs-la` for the root solve. The live
`fs-la` factorization entry points are dynamically sized and do not provide the
required finite/conditioning refusal, while taking that dependency would pull
the full GEMM/executor/allocation runtime closure into this bounded L3 crate.
The private six-by-six checked solve is therefore the narrower layer-preserving
boundary.

## Conformance tests

The inline test module in `src/lib.rs` contains nineteen focused checks:

1. mass, inertia, and centre-of-mass-reference refusals;
2. canonical quaternion double-cover selection and a known 180-degree rotation;
3. right-handed 90-degree z rotation plus nontrivial quaternion composition;
4. scaled normalization for finite huge/subnormal quaternion, axis, and
   rotation-vector inputs;
5. analytic constant-gravity translation plus energy conservation;
6. constant force/torque updates in their declared frames;
7. torque-free spherical-body energy and world angular-momentum conservation;
8. axisymmetric Euler-top phase sign and refinement toward its analytic body
   angular-velocity solution;
9. roundoff-limited energy defect for an asymmetric free body at two timesteps;
10. cancellation before a whole step, preserving the last completed state;
11. refusal of overflowing diagnostics, including an already-cancelled advance,
    without a false successful diagnostic receipt.
12. body/world point-arm transforms and `v_com + omega cross r` kinematics;
13. a rotated-arm impulse oracle for both linear and angular momentum plus its
    midpoint kinetic-work identity;
14. force-duration conversion matching the same declared impulse event;
15. a closed-form offset directional-effective-mass oracle and zero-direction
    refusal;
16. equal-and-opposite two-body impulse balance and paired work ledger;
17. non-finite and overflowing event-input refusal without a partial state;
18. bit-exact canonical-quaternion replay without renormalization drift;
19. bounded-energy symplectic kick-drift stepping under a recomputed
    harmonic force.

`tests/energy.rs` checks asymmetric torque-free energy and world angular
momentum over 50,000 steps, work balance under changing prescribed body torque
and world force, second-order attitude and momentum convergence against the
closed-form symmetric Euler top, time reversal, and exact large-angle
single-axis spin. These test the actual rigid-body owner and its real math and
geometric-algebra dependencies.

The `articulated` module additionally checks spatial-inertia matrix/direct
momentum agreement, physical-inertia and derived-overflow refusal, two-link
transform order, a closed-form gravity pendulum, prescribed base acceleration,
single-body ABA, coupled RNEA/ABA round trip, linear-storage metadata, and
position/speed/effort limit refusal.

The free-floating battery additionally checks a zero-DoF rigid body under
uniform gravity, torque-free instantaneous force/energy balance, a closed-form
external-wrench acceleration, analytic base/joint reaction coupling, common
world-rotation equivariance, a free-ABA to prescribed-base-RNEA round trip with
zero root wrench, redundant-root, non-symmetric, and ill-conditioned-system
refusal, an off-diagonal 6x6 solve oracle, the explicit spatial-to-Cartesian
origin-acceleration conversion, and linear-storage metadata. The Unitree
G1-derived 29-DoF catalog additionally proves the stronger free-fall invariant:
zero internal joint acceleration and the same world gravity expressed in every
link's body frame.

The `robot_models` module additionally checks both catalogs' link/DoF counts
and stable source order, bilateral G1 neutral-origin symmetry, an independent
iiwa neutral endpoint, admitted hard limits and physical inertias, zero dense
generalized-matrix entries, zero-input ABA, deterministic rebuilds, and retained
provenance/omission records.

The `goddard` module additionally checks its torque-free quaternions, exact RPM
conversion, Claim 2 signed ratio margin, auxiliary firing order, gyroscope
omission probe, and non-finite/out-of-domain refusal.

The `planar_drive` and `roomba` modules additionally check analytic straight,
spin, and constant-curvature motion; invalid geometry/time refusal; optical
surface-absence redirect; Claim 1 subsystem inversion; room and low-solid
non-penetration; bounded collider admission; and bit-identical replay.

These are local G0-style checks. They do not constitute contact or constrained
dynamics validation, full robot-model validation, performance evidence, G4
fault injection, or G5 cross-ISA evidence.

## No-claim boundaries

- No contacts, impacts, friction, complementarity, penetration handling, or
  static-friction capacity are represented.
- No holonomic loop constraints, nonholonomic rolling constraints,
  generalized-alpha lane, general Newton/Krylov API, or constraint
  impulse/defect receipt is present. The articulated tree supports only fixed
  and scalar revolute/prismatic/helical joints.
- The articulated lane provides a full 6x6 spatial inertia and an articulated
  tree, but no geometry-derived inertia assembly, flexible body, closed-loop
  graph, or Euler-disc-specific rule. The separate point/impulse API admits
  fixed full rotational inertia about its centre of mass, without an arbitrary
  offset-reference spatial-inertia state.
- Catalog entries are transcriptions and explicit reductions, not complete or
  manufacturer-validated digital twins. The G1 entry omits all arm DoFs and
  links, fixed body/sensor attachments, contact geometry, and their mass rather
  than lumping omitted mass into retained links. The iiwa entry omits the world
  link, fixed massless flange, damping, soft safety limits, and payload/tool.
  Both catalog roots are fixed relative to their declared base frame and can be
  passed either to prescribed-base dynamics or the free-flight base solver; the
  catalog builders themselves do not choose a dynamics boundary.
- No catalog entry loads meshes, supplies collision/contact geometry, models
  actuators/gearing/transmissions, certifies hardware safety envelopes, parses
  URDF/Xacro at runtime, or claims the omitted upstream structures have no
  physical effect.
- No collision detection, signed gap, support mapping, common-point proof,
  contact selection, impact/restitution law, complementarity solve, friction
  cone, or no-slip constraint is implemented. An equal-and-opposite impulse is
  algebraic action/reaction only; it does not by itself establish angular-
  momentum conservation about a shared contact point or physical admissibility.
- The Roomba room boundary is a display-only kinematic projection and reports
  only which declared boundary was projected plus its geometric normal. It
  supplies no contact time, impulse, normal force, restitution, friction,
  traction, wheel slip, motor torque, battery load, dust pickup, path coverage,
  localization accuracy, or hardware-validation claim. Its expanding spiral
  and randomized turn duration are contextual museum motion, not printed
  performance promises of US 6,594,844.
- Free-floating dynamics means unconstrained free flight only. It supplies no
  ground, support, contact, impact, friction, complementarity, buoyancy,
  aerodynamic, or controller force implicitly. A gravity-only robot therefore
  falls freely; it does not stand, balance, or behave as a ground-supported
  pendulum without explicit external/contact forces from another owner.
- The midpoint update is not claimed symplectic, variational, adaptive,
  adjoint-capable, or physically validated. Its torque-free energy and spatial
  momentum invariants hold to solve tolerance and roundoff for admitted fixed
  full inertia through its principal-frame map; prescribed torque obeys the
  stated midpoint work balance.
  Energy conservation alone does not bound phase/orientation error: timestep
  refinement remains necessary. The separate kick-drift lane retains its own
  first-order rotational approximation and does not inherit these invariants.
- Diagnostics describe the simulated smooth state; they do not mint a
  certificate, authority, or release-level conservation claim.
- The source-bounded Goddard apparatus lane claims no absolute patent
  dimensions, mass properties, burn rate, force, thrust, Mach number,
  trajectory, aerodynamic stability, gyroscope torque capacity, or empirical
  validation. Its ideal zero-world-rate camera support is a kinematic teaching
  limit, not a prediction that the printed mechanism rejects arbitrary loads.
- `goddard::step_goddard_rocket` uses an adjacent liquid-propellant/de Laval
  interpretation. No consumer may present those outputs as a claim, disclosed
  embodiment, or source-derived parameter of US 1,102,653. It is an ideal
  calorically-perfect nozzle calculation under declared inputs, not a
  propellant model: no chemistry, temperature-dependent heat capacities,
  real-gas behavior, boundary-layer loss, shock/separation, tank depletion,
  vehicle dynamics or empirical performance accuracy is claimed. The general
  source-backed thermochemical nozzle belongs to `fs-material`, not this
  museum teaching surface.
- The source-bounded `otis` lane composes twelve generic scalar joints for
  platform D, safety bar F, levers E and pawls f, winding drum H, shaft I,
  power drum N, shipper S, brake linkage X/Y/Z, and counterpoise R. It enforces
  the printed belt/stop/hook-lock topology but claims no historical mass,
  speed, spring rate, force, stopping distance, engagement time, power, or
  unprinted dimensions. Its lower-limit threshold is a declared normalized
  display boundary.
- The source-bounded Salisbury lane claims no historical pulley or link
  dimensions, cable material or diameter, mass, inertia, damping, motor rating,
  friction coefficient, contact modulus, grasp force, speed, stability margin,
  or force-closure guarantee. Its normalized display pose belongs to the
  browser presentation layer and is not returned by `fs-mbd`.

## Stationary-base generalized mass actions

`articulated::stationary_base_joint_momentum` evaluates `M(q) qdot` through the
existing zero-bias RNEA recursion. `stationary_base_joint_velocity` applies the
inverse mass through the existing ABA recursion, retaining linear storage and
physical-inertia checks. Neither builds a dense generalized matrix. Momentum
has units N s or N m s for prismatic or revolute coordinates; it is not an
actuator force and does not use force limits. Actual input/output positions and
rates retain their native joint-limit checks.

These functions assume a stationary prescribed base; a moving base adds an
affine momentum term and requires a separate boundary. Their immediate consumer
is source-heated articulated dynamics in fs-mbd-material, where thermal inertia
updates preserve joint momentum and account parameter work explicitly.
`cargo test -p fs-mbd --test articulated_mass` uses an independently known
two-joint inertia matrix and checks the momentum/force-limit distinction.

## Free-floating articulated impulse response

`articulated::free_floating_momentum` sums every link's kinetic energy and
world-frame linear plus spin/orbital angular momentum using full offset spatial
inertia. `apply_free_floating_impulses` solves the coupled six-base-plus-joint
inverse-mass action through the existing O(n) articulated-body algorithm with
velocity and gravity bias zero. Per-link body impulses act about link origins;
internal generalized joint impulses recoil the free base without creating net
world momentum. Poses and configuration remain unchanged. Actual output speed
limits are checked without clipping, and impulse values are never compared to
actuator-force limits. The atomic result reports both endpoint momenta and
energies, total external world impulse, midpoint impulse work, and kinetic-work
residual. This owner supplies inertial impulse response; contact selection,
restitution, friction, CCD, and time integration remain separate responsibilities.

G0/G1/G3 `tests/articulated_mass.rs` checks the independent offset-COM body
solution, full three-inertia coupled response, internal reaction momentum,
world-frame covariance, impulse-work closure, and atomic admission.

The impulse receipt also exposes actual linear/angular momentum increments and
`kinetic_energy_change`, evaluated from the returned link velocities at the
fixed configuration. Kinetic change uses the polarized COM energy form; spin,
COM-arm and world-origin orbital momentum increments are accumulated separately
with the existing double-double scalar owner. This retains small physical work
and recoil alongside a large unchanged translation or orbital momentum. The
reported increment is not the requested load or an inverse-mass probe: a load
too small to change the stored endpoint remains visible in the work residual.
Free and prescribed impulse owners now use the shared original spatial-inertia
kernel before rounding differences, COM velocity cross terms or link sums.
Independent impulse work retains the mean-rate product layers through its
comparison with storage. The public `retained_kinetic_energy_change_j` and
`retained_impulse_work_j` accessors expose the retained scalar values; supported
events also expose `retained_support_work_j`. Their f64 fields remain diagnostic
roundings. Single-plane and shared normal/frictional/stateful contact mechanics
retain the native product layers across systems and source work. Native velocity
updates, endpoint speed limits and constitutive tolerances remain unchanged.
Focused regressions cover boosted contact loss, offset-world angular impulse,
internal joint work and recoil, an unrepresentable requested velocity increment,
and a partially rounded velocity jump whose independent work defect remains
nonzero even when the printed kinetic change and work agree.

## Continuous free-floating articulated trajectories

`articulated::motion::step_free_floating` advances an actual base pose/twist and
all joint coordinates/rates under held internal efforts, held link-frame
wrenches and uniform world gravity. `FreeFloatingMotionState` composes the
canonical free-base state with the existing compact joint coordinates. Its
accepted output is directly reusable as the next input; the function never
mutates the caller's state or publishes a partially advanced result.

The second-order explicit Lie-midpoint lane integrates total world linear
momentum, angular momentum about the system COM, COM position, base orientation,
and joint position/rate. Native free-base ABA supplies joint accelerations.
O(n) shape kinematics assembles the composite locked inertia at the COM and the
momentum due to relative joint motion. The existing fixed 6x6 physical-inertia
solve recovers base velocity at the supplied joint rates. This is a coordinate
inversion of the momentum map, not a velocity clamp or a solve with an
articulated inertia that permits extra joint response. Reaction wheels and
extending limbs therefore recoil their base without creating total momentum.
No dense generalized mass matrix is constructed; work and storage remain O(n).

The midpoint COM is `Xh = X0 + dt P0/(2 M)`, while the accepted COM uses
`X1 = X0 + dt (P0+P1)/(2 M)` and `P1 = P0 + dt Fh`. Thus the discrete orbital
angular change is `dt Xh cross Fh`. Advancing central spin by the midpoint COM
moment preserves the full world-origin angular impulse of the same stage
wrenches, up to numerical roundoff. Lever arms and inertia are assembled near
the COM before world translation, avoiding subtraction of large orbital
momenta. Body forces rotate with each link at the force stage; gravity acts at
the total COM and cannot generate a spurious central moment.

The receipt independently evaluates actual endpoint physical momenta and kinetic
energy, reports the external world impulse, and separates actuator work from
external work. Uniform gravity uses exact COM displacement work; other work is
midpoint quadrature. Its kinetic/work residual is integration error, not physical
dissipation. This lane is neither energy-exact nor symplectic, and does not
supply contact, loop constraints, adaptivity or a controller. Step refinement
remains necessary for fast internal dynamics.

The lane admits at most 4096 links so each indivisible native tree pass is
bounded. Cancellation is polled between those passes and every 64 links in new
momentum passes, including a final poll before publication. The native fixed-root,
physical-inertia, finite-value, and joint position/rate/effort bounds remain in
force at the initial, midpoint and accepted states; limits are never clipped.
`tests/articulated_motion.rs` checks analytic reaction-wheel and telescoping-mass
trajectories, full off-center force/angular impulse, branched mixed-joint free
momentum and COM drift, second-order convergence to the exact symmetric-top
solution, frame covariance, distant admitted origins and atomic refusal/cancellation.

### Conservative spatial gravity on original articulated links

`articulated::motion::gravity::PreparedArticulatedGravity` binds the original
tree and prepares each link's complete central inertia once through the native
rigid mass-property admission. Uniform, central and quadrupole fields sample
the actual COM poses at the existing initial and midpoint motion stages. COM
forces and tidal torques shift to the original link origins for ABA, while
their physical COM resultants and power drive the same O(links) momentum owner.
The separate held uniform gravity must be zero. Endpoint field admission and
per-link potential changes use the already computed native endpoint kinematics;
any refusal returns no candidate. Gravity work remains part of the original
external-work receipt, with a separate potential account and measured mechanical
work defect; this does not assert exact discrete energy conservation or source
recoil. `tests/articulated_gravity.rs` checks uniform/free-fall and translating
frame equivalence, differential tidal joint acceleration, full-tensor torque and
COM-wrench potential derivatives, offset-COM Kepler refinement, independent
nongravitational work, and late field-domain refusal.

### State-dependent internal actuator effort

`articulated::motion::step_free_floating_with_effort` evaluates a caller's joint
effort law at the actual initial and reconstructed midpoint states, at elapsed
times zero and `dt/2`. Each output slice starts from the admitted held efforts,
so motor torque-speed curves, torsional springs/dampers and smooth controllers
can add to or replace them. The original ABA and locked-inertia momentum map
retain the full carrier reaction; internal motor torque never becomes external
angular impulse. The exact accepted midpoint effort and relative joint rate
supply actuator work. Revolute/helical effort is in N m; prismatic effort is in N.

`spatial_loads::step_free_floating_spatial_with_effort` composes the same law
with original moving-point spatial loads and finite axial connections at those
same stages; `spatial_loads::bushing::step_free_floating_bushings_with_effort`
also retains anisotropic finite-rotation bushing forces, storage and damping.
Held APIs retain their original arithmetic and borrow their effort
slice; driven APIs reuse one O(joints) buffer. Finite-value and native effort
limits are checked at both stages, and cancellation brackets each callback.
Inputs remain immutable on every refusal; caller-owned callback side effects
are not rolled back. Second order applies to smooth resolved load laws. This
explicit lane does not imply stiff stability, motor electrical-state coupling,
hybrid controller switching or a dissipation/thermal model for callback effort.

### Stage-dependent external forcing and internal actuation

`articulated::motion::step_free_floating_with_forcing` samples external link
wrenches and internal joint efforts together at elapsed times zero and `dt/2`,
using the actual reconstructed stage states. Both output slices reset from the
held loads before each call. External wrenches remain about link origins in
their own body frames; uniform gravity retains its separate original owner.
The same accepted samples drive native ABA, momentum reconstruction and the
original external impulse/work and internal actuator-work accounts. Internal
efforts retain finite carrier recoil and never enter external momentum supply.
One wrench and one effort buffer are reused; existing held/effort-only paths
retain their arithmetic. Callback errors, nonfinite loads, effort limits and
cancellation refuse atomically. Smooth histories and feedback must be resolved
by this explicit step; load jumps require step splitting, impulses have their
separate owner, and callback work is not automatically assigned to heat.
Focused tests in `articulated_motion` cover ramped thrust, analytic drag decay,
combined motor/external torque, held-path equivalence and midpoint refusals.

## Live free-floating inertia reconfiguration

`articulated::motion::reconfigure::transfer_inertia` transfers an actual moving
tree between full centroidal link tensors at unchanged configuration. It
preserves all six base and every joint conjugate momentum:
`M_next(q) u_after = M_previous(q) u_before`. The original articulated-body
inverse-mass recursion applies the individual old link momenta to the new
operator, with O(links) work and storage and the existing fixed-size root solve.
There is no assembled generalized mass matrix or fictitious external impulse.

Link identities, topology, transforms, joints and limits remain identical.
Per-link mass and COM must be bit-identical. The caller supplies the
closed-inventory interpretation; this reduced parameter update does not infer
internal shape velocity/momentum, moving attachments, mass flow, contact,
closed-loop constraints or fluid storage.

The result independently reconstructs all 6+n momenta through a zero-velocity
RNEA mass action and world linear/spin-plus-orbital momentum through native link
kinematics. Caller-declared componentwise physical tolerances govern admission.
The inverse-mass solve uses a common translating frame to avoid contamination
by large bulk motion. Signed parameter work uses compensated COM-velocity
differences and the full centroidal-tensor identity
`delta Krot = mean(omega) dot I0 delta(omega) + omega1 dot (I1-I0) omega1 / 2`.
Every original link enters the shared kinetic accumulator before rounding,
retaining small work beside unchanged bulk translation or spin. The actual
accepted endpoint governs the absolute work budget, for positive and negative
work alike. No heat partition follows.

Identical models preserve state bits after admission. Invalid resulting rates,
changed inventory/geometry, work exhaustion, numerical failure or cancellation
publish no state. The native 4096-link cap bounds indivisible tree passes.

Three focused native regressions cover an analytic carrier/rotor transfer and
continued flight, a mixed-joint offset tree with full tensors and independent
link-level momentum reconstruction, frame/Galilean covariance, reverse transfer,
large common translation, and atomic identity/limit/work/cancellation behavior.

## Rigid payload capture on a moving articulated tree

`articulated::motion::inventory::capture_payload` attaches an independent native
spatial-inertia body to any existing link at its actual supplied pose. The link
frame, topology, joint coordinates and base pose remain unchanged. Loaded mass,
COM and full rotated centroidal inertia use a reduced-mass parallel-axis sum;
a mass increment below representable resolution refuses. The original native
articulated inverse-mass recursion solves all base and joint recoil from the
combined incoming momenta. No dense generalized mass matrix is assembled.

Optional existing constant affine gear/lock rows react with the loaded inertia
before final speed-limit admission. Their initial position and velocity must
already satisfy the declared rows. The bounded prepared row-mobility owner is
shared with constrained inertia reconfiguration; there is no independent gear
solver or intermediate free-candidate speed rejection.

The receipt retains the loaded native model/state, the actual attached payload
state and an opaque attachment binding the original retained inertia and exact
loaded target signature. Independent world linear/angular momentum sums and
zero-velocity RNEA reconstruct every base/joint conjugate momentum, including
the incoming payload pulled back to the target's virtual link. Only the
reported internal gear reaction changes total projected joint momentum.

Removed relative energy is a sum of nonnegative old-body velocity-change
quadratics. The inverse-mass solve and separate endpoint energy audit use one
common translating inertial frame, preserving small capture signals under a
large uniform velocity. Actual finite gear midpoint work is reported and capped
separately. No rounding defect is assigned to heat. All physical momentum,
energy-loss, reaction-work and balance tolerances must pass before publication;
cancellation and native inertia/state/limit refusals leave inputs untouched.

Three focused native regressions cover analytic carrier/rotor capture and
continued native motion, a 3.8475 J loss under a billion-scale common velocity
with a 4 J allowance, mixed revolute/prismatic/helical joints with offset COMs
and full rotated tensors, independent spatial/screw and parallel-axis oracles,
scene/boost covariance, locked-joint recoil with an otherwise inadmissible free
candidate, and atomic inventory/budget/constraint/cancellation refusals.

This is a declared instantaneous rigid attachment. Contact discovery, finite
grasping duration, fracture, mass flux, internal shape momentum and the thermal
destination of removed energy require their own physical models.

## Articulated payload release with internal separation recoil

`articulated::motion::inventory::release::release_payload` consumes the capture
attachment and current loaded tree. It restores the target's exact retained
spatial inertia and reconstructs the independent payload at the current link
pose with its actual inherited rigid velocity field. The target's name, parent,
reference transform, joint and loaded inertia must still match the attachment.
Other links keep their current inertias. Nested captures on one link release in
reverse inventory order; repeated release against restored inventory refuses.

An optional `SeparationImpulse` specifies one link-frame application point,
force impulse on the payload and free angular impulse. Its point-force moment
and free couple form one complete wrench; the retained link receives the exact
opposite wrench. The original articulated inverse-mass owner supplies coupled
tree recoil, and the original six-component spatial-inertia solve supplies
payload recoil. Optional ideal gear/lock rows react with the retained inertia
before endpoint speed admission. Zero impulse still admits the state and row
rank, then preserves the retained tree's state bits.

Both surviving native states, the actual impulse pair and separate separation
and gear work are returned. Separation work pairs the full impulse with actual
midpoint endpoint velocities, including force/couple and constraint recoil
cross terms. Common-frame energy, world momentum and independent projected
total momentum audits use the actual outputs. The released payload is pulled
back to the target's virtual link solely for the fixed-event projected audit;
internal separation impulses cancel and only gear reactions change total
projected joint momentum. Absolute separation-work and finite gear-work
allowances are independent, and numerical residuals remain explicit.

Capture and release retain `kinetic_energy_change_j` from every original
retained link and physical payload through the shared compensated spatial
kernel. Absolute common-frame energies remain diagnostics. Repeated-target
inventory prefixes telescope at one accepted endpoint, without sequential
dynamics. Supported entries retain the same calculation in laboratory and
accounting frames.

`inventory_decomposition_j` separately reports the signed finite difference
between assembled and constituent storage. Original spatial coefficients are
differenced before the velocity quadratic; the native payload adjoint and a
polarized correction retain its actual rounded endpoint twist. Capture and
release apply opposite decomposition signs. This gated numerical storage defect
is never capture heat, release work or support work. The supported receipt also
retains `accounting_inventory_decomposition_j`. Added accounting is linear in
links plus payloads, with fixed stack storage per local decomposition.
Three executed native kernel groups cover small signed exchange beside dominant
spin, unresolved inertia coefficients, full tensors, offset COMs and payload
transform rounding. Four integration groups specify free/locked, batch and
driven-support lifecycles and assembly-defect refusal; those complete native
integration tests remain unexecuted in this environment.

Three focused native regressions verify zero-impulse release after native
motion and continued motion of both outputs; a combined point-force/free-couple
case with 0.99 J work under a billion-scale uniform velocity; locked-joint recoil
with 16/75 J work and an inadmissible free intermediate; independent full-tensor
principal-axis payload recoil and shared-point work; physical momentum/energy
and frame/boost covariance; stale inventory, work limits and late cancellation.

The caller supplies the actual release event, attachment geometry and work
source. This instantaneous model does not infer fracture, finite-duration
separation, subsequent contact, mass flow or thermal conversion.

## Continuous prescribed-inertia motion

`articulated::motion::reconfigure::step::advance` and
`advance_with_effort` consume full centroidal-tensor samples at the start,
midpoint and end of a window. They preserve every base/joint momentum while
transferring to the midpoint model, call the original native free-floating
smooth owner, and transfer at the actual smooth endpoint to the final model.
The original articulated recursion and explicit-midpoint dynamics remain the
single numerical owners; no dense generalized mass solve is introduced.

The same per-link mass, COM, topology, attachments, joints and limits are
required throughout. This is prescribed reduced inertia with no additional
intrinsic shape momentum, mass flux or moving-COM degree of freedom. Smooth
tensor history, midpoint samples accurate to O(duration²), smooth resolved
loads and adequate mode resolution give second-order refinement. Discontinuous
changes use the separate transfer event. Symmetric placement of the transfers
does not make the original explicit-midpoint step reversible, variational or
energy-exact.

Both transfer receipts and the complete native motion receipt are retained.
Signed parameter work remains separate from actuator, external and internal
body work. The complete-window budget caps the sum of absolute parameter work:
the final transfer receives only the unused allowance, so opposite-signed
changes cannot cancel out of admission. Actual endpoint energy and world
momentum residuals independently audit the entire step. The public
`kinetic_energy_change_j` is independently reconstructed from accepted endpoint
rates and full tensors. Storage is compared with the original parameter,
actuator and external/internal-body work channels before rounding, rather than
subtracting absolute kinetic totals or defining a source by residual. The fixed
per-link mass/COM and prescribed-inertia model restrictions remain unchanged.

Stage-dependent effort callbacks retain the native elapsed times zero and
half-duration. The first callback sees the virtual transferred state at the
midpoint inertia, and the second sees the native predicted midpoint state.
No callback or elapsed time is introduced by a parameter transfer. Invalid
models/duration, either transfer or native refusal, exhausted work or
cancellation publish no state; external callback side effects cannot be rolled
back. Hidden actuator storage, heat and dissipation need their own models.

Three focused native regressions cover driven analytic carrier/rotor rates and
momenta, second-order angle and parameter-work convergence, exact unchanged-model
equivalence, a zero-net-work cycle with 28/15 J gross work, pre-callback admission,
midpoint cancellation, and state-dependent damping with independently integrated
mechanical and actuator work. Refinement uses 8, 16 and 32 steps over one second.

## Gear-constrained live inertia changes

`articulated::motion::reconfigure::constraints::transfer_inertia_constrained`
accepts constant affine joint relations `A q = b`, including ideal gear ratios
and locked coordinates. Per-row coordinate units, coefficients, position and
rate tolerances are explicit. Initial configuration and tangent rates must
already be admissible. The original fixed-mass/COM, unchanged-topology tensor
compatibility rules remain in force.

The coupled event solves `M_next u_next = M_previous u_previous + [0; A^T lambda]`
with final `A qdot = 0` to the declared tolerance. Original native ABA mass
actions form a diagonally normalized mobility for at most 32 independent rows;
the existing prepared LU owner solves it and retained-factor unit-RHS solves
measure numerical reciprocal conditioning. There is no assembled generalized
mass inverse or regularization of dependent constraints. With fixed row count,
mass actions and response storage scale linearly in link count.

Only the complete constrained endpoint receives native speed and work admission.
An inadmissible unconstrained intermediate therefore cannot reject an admissible
coupled result. Independent generalized-momentum audits include the actual
reaction covector, while native world linear/angular momentum must remain
conserved. Conditioned kinetic change, parameter work and actual midpoint
reaction work are separate. Both work channels have independent absolute caps;
finite tangency residual work is retained and no impact heat is inferred.

Empty rows exactly delegate to the original transfer. Identical models preserve
state bits after constraint and mobility admission. Invalid geometry/rates,
rank or conditioning, exhausted work or cancellation publish no state.
The G0/G1/G3-focused tests use analytic geared rotors, an offset full-tensor
mixed-joint tree with independent momentum reconstruction, frame/boost and
signed-row-scaling transformations, reverse transfer, and a locked-joint case
whose free intermediate violates limits while its coupled endpoint is valid.
This slice covers constant ideal joint relations; moving ratios, backlash,
nonlinear point/weld closures, rolling and unilateral contact need their own laws.

## Continuous gear-constrained inertia motion

`articulated::motion::reconfigure::constraints::step::{advance, advance_with_effort}`
surround the original native free-floating step with constrained start-to-midpoint
and midpoint-to-end tensor transfers. At both actual native force stages, the
original free acceleration and the mobility at that stage's configuration supply
joint reactions enforcing `A qddot = 0`. Constant affine gear/lock relations
therefore remain compatible with the native midpoint position/rate updates.
Initial, stage and final row residuals are admitted; no position repair is added.

Caller effort is evaluated at native elapsed times zero and half-duration.
The first callback sees the virtual constrained rates under midpoint inertia.
One original ABA evaluates each caller effort; a native inverse-mass response
adds the ideal reaction directly to its joint acceleration. There is no second
dynamics pass. Motor effort limits admit the caller alone, while real state
limits still apply to the resulting constrained motion. The native receipt
now separates `actuator_work` from signed `constraint_work`, both using the
same actual midpoint rates. Passive reaction work never becomes actuation or
damping loss.

Separate complete-window caps count absolute parameter and reaction work.
The last transfer receives only each remaining allowance, including reaction
work already supplied during smooth motion. Opposite signs cannot hide work.
Both transfer receipts, both stage reactions and the native motion result are
retained, with independently calculated final energy and world-momentum residuals.
The same compensated endpoint-storage owner supplies `kinetic_energy_change_j`
and the full-window residual. Parameter, transfer-reaction, smooth-reaction,
caller-actuator and external/internal-body work remain distinct inputs to that
comparison.
A failed stage/helper/native step or cancellation returns no state; typed helper
errors survive the callback adapter. Callback-owned side effects cannot roll back.

Second-order refinement requires smooth resolved loads and prescribed tensors,
an O(duration²)-accurate midpoint sample and adequate mode resolution. This
does not turn native explicit midpoint into an energy-exact or reversible method.
The focused G1 tests compare driven nonuniform inertia against scalar reduced
gear equations and analytic reaction torques, then state-dependent damping
against closed-form motion and separately integrated work. Both use 8/16/32
steps over one second; empty-row equivalence, gross work and midpoint cancellation
are also exercised. A locked finite slider independently checks 60 N passive
reaction with a 0.1 N motor rating, complete carrier recoil and external work;
a finite tangency defect checks separate signed reaction work and its budget.
This remains a constant ideal-joint model, excluding
time-dependent ratios, constraint engagement, backlash and general contact.

## Full-inertia rigid RATTLE projection

`rigid_rattle::PointPlaneRattle` implements the canonical
`fs_time::se3::RattleProjection` seam for one fixed world plane acting on one
material point of the original COM-referenced rigid body. The scalar constraint
is `g = n·(x + R r - p)=0`; in body coordinates its Jacobian is
`J=[r×n_b,n_b]`. Both configuration and tangent corrections use the actual
kinetic metric `diag(I_body,mI)`: `delta = M^-1 J^T lambda`. Full products
of inertia therefore affect how an offset-point correction divides between
rotation and translation. There is no diagonal-inertia shortcut, reference-point
clamp, penalty stiffness or velocity clipping.

Position projection uses bounded Newton corrections on SE(3) and independently
rechecks the actual material-point gap after every accepted correction. Velocity
projection is the exact scalar mass-metric tangent projection at the supplied
pose and rechecks the actual point normal speed. The implementation exposes the
physical scalar mobility `J M^-1 J^T` and a receipt containing removed position
and velocity violations plus final residuals. Invalid geometry, a zero normal,
nonfinite response or exhausted nonlinear iterations refuses without publishing
a partial candidate.

This projection interface is a reusable constrained-state repair. It does not infer unilateral activation, restitution,
friction, support work, continuous contact force, multiple simultaneous
constraints or event timing. The caller selects when a bilateral point/plane
constraint is physically active and composes this projection with an
`fs-time` variational/SE(3) step. The two focused tests verify a nonprincipal
full tensor against the independent momentum-covector condition
`M delta = J^T lambda`, show that products of inertia create cross-axis angular
correction, verify nonlinear gap closure and tangent closure, and cover explicit
invalid-geometry refusal.

## Constrained rigid point/plane motion

`rigid_rattle::PointPlaneRattle::advance` evolves an already constrained native
rigid state using an initial physical normal impulse, the original full-inertia
free step, and a normal impulse at the actual endpoint material point. The
returned pose is the free integrator's endpoint; it is never repaired. Under
compatible external forcing, the plane's translation and rotation symmetries
retain tangent linear momentum and plane-normal total spin-plus-orbital angular
momentum. Both signed impulse works remain visible in the returned receipt.

The scalar position solve starts from zero impulse, uses the original normal
mobility and accepts only residual-decreasing corrections with bounded
backtracking. Explicit impulse, angular-travel and iteration caps admit a local
smooth branch; failure requires physical timestep refinement or revised bounds.
Initial gap and normal speed must already meet the projector's tolerances.
Cancellation or a numerical refusal publishes no state, and the new step
allocates no heap storage. The original world-fixed force/body-fixed torque
convention is retained.

`tests/rigid_rattle_step.rs` checks exact COM support under gravity, full-tensor
offset-point momentum and impulse-work balance, reversal, second-order
refinement, rigid-frame/tangent-boost covariance and bounded refusals. This is
a fixed-plane bilateral law, including tensile reactions. It does not claim
exact energy conservation, symplecticity, global root uniqueness, unilateral
contact, friction or moving supports.

### Prescribed translating supports

`PointPlaneRattle::advance_translating_plane` composes the same native impulse,
free step and endpoint impulse with prescribed normal translation. A
`step::PlaneTranslation` supplies initial normal speed and held normal
acceleration; the plane point in the original constraint is the initial point.
The position solve targets `p+n*(v*dt+a*dt²/2)` and the endpoint speed constraint
targets `v+a*dt`. Full-inertia, off-center coupling and the original local
impulse/angular-travel/backtracking controls remain shared with the fixed law.

The receipt returns the advanced plane and speed for continuation, both native
impulse events, their reciprocal support impulse, and signed boundary work.
Each reaction is paired with the prescribed support speed at its own event.
The independently computed relative-frame work defect stays numerical; it is
not heat or an energy correction. Prescribed supports supply external momentum
and power; no finite support mass or inertia is inferred.

Endpoint velocity and displacement are accumulated in double-double arithmetic
before rounding. Their representability residuals and combined endpoint errors
must fit the original physical tolerances. A small acceleration hidden beside a
large carrier speed therefore refuses when it cannot be represented. No heap
storage or partial state is published. The stationary motion recovers the
original native result exactly. Plane rotation, friction, unilateral contact
and finite-support recoil remain outside this entry.

`tests/rigid_translating_plane.rs` checks analytical accelerating COM support
and signed work, oblique normals, full-tensor offset-point boost covariance,
reversal and second-order refinement, continuation, momentum symmetries,
stationary identity, unresolved prescribed acceleration and all cancellation
checkpoints.

### Prescribed rigidly rotating supports

`PointPlaneRattle::advance_rotating_plane` extends the same native solve to a
plane rotating at a held world angular velocity about a translating pivot.
`step::PlaneRigidMotion` declares the pivot, its held world velocity and the
held world angular velocity. The pivot need not lie on the plane. The endpoint
normal and plane point follow the original native quaternion rotation; the
velocity constraint uses the actual source field `U + omega cross (x-c(t))`
at the material point. The initial and final physical impulses use their own
endpoint normals, with no position repair or substituted free-body dynamics.

The receipt returns the endpoint plane and motion for continuation, reciprocal
linear impulse, reciprocal angular impulse about the fixed INITIAL pivot,
signed support work, and a separately paired relative-frame work defect.
Prescribed motion supplies external momentum and power; no finite support
inertia is inferred. Tangential surface motion contributes no frictionless
normal work. Constraint work and source work are not relabeled as heat.

Double-double geometry retains prescribed motion until physical admission.
Actual material-point offsets use `(COM-origin)+arm_world`, and relative point
velocities retain the separate native COM and spin components. A remote world
origin or large common boost therefore cannot erase a short lever arm, local
normal speed or torque. Both represented and pre-rounding prescribed endpoint
gap/speed must pass the original tolerances; an unrepresentable moving pivot
refuses rather than changing the source at continuation. Static/translating
entry arithmetic remains unchanged, and stationary rigid motion recovers the
original step within this admitted physical geometry domain.

The original impulse, iteration and backtracking caps remain authoritative;
both body and prescribed support travel must fit the local angular cap. The
step owns no heap storage and cancellation or late refusal returns no state.
This remains a smooth bilateral local branch, including tensile reactions,
without unilateral contact, friction, a finite-inertia support, exact energy
conservation or a global nonlinear-root certificate. The G0/G1/G3/G4 cases in
`tests/rigid_rotating_plane.rs` cover an independently derived discrete
point-mass solution, the continuous centrifugal trajectory and second-order
refinement, reversal, full-inertia frame/boost covariance, P/L/source-work
closure, stationary identity, remote-origin and hidden-spin failures, and
every cancellation checkpoint.

### Plane constraints carried by a finite rigid body

`rigid_rattle::finite_plane::FinitePointPlaneConstraint` constrains a material
point of one native rigid body to an infinite plane carried by a second finite
body. Both evolve with their actual full inertia, original native free motion
and independent applied wrenches. The normal reaction acts on the carrier at
the current projected sliding point. Its equal-and-opposite impulse therefore
retains the carrier's rotational Jacobian and reciprocal torque; a fixed plane
anchor is not substituted for the physical reaction point.

The scalar impulse/free-step/impulse engine is shared with the original rigid
distance pair. Its mobility arithmetic, native reactions, local backtracking,
angular/impulse limits and cancellation order retain the original distance
path. Initial position and velocity constraints must hold; neither pose is
repaired. Every accepted receipt returns both bodies, both native pair events,
both free steps, gap/speed residuals and numerical constraint work evaluated
from relative COM/spin components. This internal work defect is not an external
support source or heat. Central reactions retain total linear and angular
momentum to the admitted geometry and floating-point errors; native event
receipts retain actual momentum increments and their defects.

Double-double relative geometry preserves local arms at remote origins and
cancels common velocity before local spin is added. The plane is represented
by its original body-normal scalar offset, prepared once; translating its
authored anchor arbitrarily far along the same plane cannot introduce a false
world-normal offset or overflow an irrelevant anchor query. The reconstructed
native reaction arm must fit the original position tolerance.

The step allocates no heap storage and publishes neither body on numerical or
late cancellation refusal. This is a smooth bilateral sliding joint, including
tensile reactions, without unilateral collision, friction, a finite plane
boundary, exact-energy conservation or a global nonlinear-root certificate.
`tests/rigid_finite_plane.rs` covers analytical unequal-mass load sharing, an
independent continuous finite-inertia sliding-bead/rotor reference with
second-order refinement and reversal, full-tensor P/L/work balance, frame/boost
covariance, tangential-anchor equivalence, hidden-arm/spin admission and every
cancellation checkpoint. The original six distance-pair regressions also pass
against the shared engine.

## Distance between independent finite rigid bodies

`rigid_rattle::distance::RigidDistanceConstraint::advance` evolves two independent
native rigid bodies with a positive fixed separation between arbitrary material
points. Both finite masses recoil through native central equal-and-opposite
impulses around the original full-inertia free steps. The initial impulses close
endpoint distance; the endpoint impulses remove relative range rate. There is
no pose repair, and internal reactions preserve total world linear and total
orbital-plus-spin angular momentum. The receipt retains both pair impulse
events, both smooth steps, actual endpoint residuals and signed constraint work.

The bounded scalar solve reuses the original summed point mobility and the
existing point/plane impulse, angular-travel and work controls. Each body's
rotation and the change of connection direction must meet the local angular
cap. The inputs must already satisfy the position and velocity constraints.
Collapsed geometry, unrepresentable arithmetic, exhausted local bounds,
nonconvergence and cancellation publish neither body. Evaluation allocates no
heap storage. The force/world and torque/body conventions and shared uniform
gravity are those of the original rigid integrator.

`tests/rigid_distance_step.rs` checks analytic unequal-mass dimer recoil and
orbital motion, differential forcing, full-tensor off-axis momentum/work balance,
reversibility, second-order phase refinement, frame and boost covariance, body
exchange and atomic refusals. This is one bilateral positive-length rod, with
tension and compression. Zero-length ball joints, general joint graphs, unilateral
cables, exact energy conservation and symplecticity remain outside this entry.

### Spatial loads during finite-pair constraint motion

`RigidDistanceConstraint::advance_spatial` and
`FinitePointPlaneConstraint::advance_spatial` compose the original scalar pair
constraint owner with `RigidBodyIntegrator::step_spatial`. World-held point
forces retain their actual rotating material arms in every trial. Body-held
and world-held free couples retain their respective frames. A once-converted
body torque would instead change the physical loading and is not substituted.

The native pair impulse/free-step/impulse solve is shared with the existing
held-wrench entry. A private statically dispatched free-step seam retains the
original `step` path and arithmetic for held wrenches, while the spatial entry
uses the original spatial-load Verlet/free-rotation owner. Its complete
`SpatialStepReceipt` is retained in each `free_steps` entry: external world
linear/angular impulse, force/couple work, and native numerical defects remain
separate from reciprocal constraint impulses and their signed work. Both public
pair receipt types have a default free-step type, preserving the original
held-wrench result layout. This introduces no alternate rigid integrator.

Prepared force moments are reused at constant cost. The orientation-independent
triangle bound from `PreparedSpatialLoads::torque_norm_bound` is computed once
per body and interval, and bounds angular travel before each free trial. It can
be conservative near a moment cancellation; a local-branch refusal requires a
smaller interval, not a frozen initial-moment approximation. There are no heap
buffers, repeated point-force reductions, or extra inertia factorizations.

Loads are held for one admitted positive interval. State-dependent springs,
drag, arbitrary constraint graphs and contact remain separate owners. Neither
exact loaded energy conservation nor continuous-trajectory enclosure is claimed.
Cancellation or late load/solve refusal publishes neither finite body.

`tests/rigid_spatial_constraints.rs` checks independent continuous loaded-rotor
and sliding-bead/finite-carrier equations, second-order refinement and reversal,
full-tensor external/reaction momentum and work accounts, frame/boost covariance,
native spatial-receipt identity, the torque bound, and all observed cancellation
checkpoints. Existing pair, finite-plane and spatial-load regressions continue
to exercise their original public paths.

## Coupled distance graphs between native rigid bodies

`rigid_rattle::distance::graph::advance_distance_graph` advances multiple
independent native bodies connected by positive `DistanceLink` rods. Shared
bodies participate in the complete `J M^-1 J^T` mobility, retaining all
translational and full-inertia angular cross responses. A coupled initial
impulse solve closes every endpoint length after the original free steps;
a coupled solve at the actual endpoint removes every range rate. No body
is fixed, no pose is repaired, and no rod is replaced by a stiff spring.

The original pair owner applies the simultaneously solved central reactions at
unchanged poses in caller link order. Both reaction inventories and every
native free step remain in the receipt. Their actual summed impulse work
includes shared-body cross terms; individual pair work allocations depend on
the reported application order. Internal reactions preserve total world
linear and orbital-plus-spin angular momentum. Finite constraint work remains
explicit and is not heat or an energy correction.

At most 32 independent scalar rows are admitted. Unit-diagonal mobility scaling
and the original fs-la Cholesky enforce the authored squared-pivot floor.
Dependent or insufficiently regular rows refuse; this is floating regularity
admission, not a condition-number certificate. Initial positions and velocities
must already satisfy the constraints. The original dimensioned impulse,
rotation, iteration and backtrack controls remain in force, with explicit
body/row, trial-evaluation and logical-payload caps. Trial storage is reused:
O(B+R) native work per trial, O(R²) mobility assembly and two O(R³) bounded
factorizations. New transaction buffers reserve fallibly; the original
bounded Cholesky allocation semantics are retained. Cancellation is checked
per body/row and before atomic publication of all states.

`tests/rigid_distance_graph.rs` checks analytic unequal-mass chain load sharing,
closed-triangle collective rotation and second-order phase, 300 full-inertia
off-axis steps with momentum/work balance, reversibility, frame/boost and
body/row permutations, pair/free/disconnected equivalence, a real 32-row
connected framework, rank/resource/local-branch refusals and cancellation
at every observed checkpoint. This entry supplies smooth regular bilateral
distance graphs; the same owner also supports spherical joints through the
mixed graph entry below. Weld graphs, redundant constraints, unilateral cables,
contact/impact coupling, exact energy conservation,
symplecticity and validated whole-path enclosures remain separate work.

### Prescribed-length internal actuators on native rigid graphs

`rigid_rattle::distance::graph::advance_driven_distance_graph` extends the same
simultaneous solver with one `DistanceDrive` per link. The prescribed length is
`L(t)=L(0)+r(0)t+at²/2`; initial range rates and endpoint lengths/rates enter the
original full shared-body mobility and native impulse/free-step/impulse path.
Mixed static/driven rods, loops, disconnected islands and unlinked bodies keep
their actual finite inertia. The original static API is unchanged. Driven
range-rate admission and work retain native COM velocity and angular-velocity
cross arm separately, including for stationary commands: a boost cannot erase
local spin by rounding a convenience point velocity first. No coordinate is
snapped or body silently converted to a support.

Central internal reactions retain total linear and orbital-plus-spin angular
momentum, while a prescribed actuator supplies or recovers signed energy.
Its work is `j_start*r_start+j_end*r_end`, with each simultaneously solved
scalar impulse paired with its own event's commanded rate. Per-link and total
work are returned separately from the actual numerical constraint-work defect.
That defect uses double-double relative material-point velocities and the
original applied vector impulses, retaining shared-body work cross terms without
subtracting large laboratory-frame kinetic energies. Per-link defect attribution
uses native application order; the total is invariant up to floating roundoff.
Neither this defect nor free-step error is physical heat.

Command preparation retains length/rate increments in double-double arithmetic.
The full prescribed quadratic, including an interior contraction minimum, must
stay above the position tolerance. Rounded continuation length/rate errors are
reported and admitted against physical tolerances; the final state must satisfy
both the original polynomial and the represented continuation constraints.
A nonzero prescribed displacement rounded to unchanged length refuses. Updated
links and rates are returned for continuation. The native body/row, impulse,
angular, iteration and trial caps remain authoritative; a bounded stack command
bank and one fallibly reserved continuation vector add O(R) work/storage. The
reported logical payload includes that vector, and cancellation remains atomic.

`tests/rigid_distance_graph_driven.rs` checks independent unequal-mass actuator
force/work balances, held external forces under a large common boost, exact
static parity, full-tensor closed-framework continuation/reversal and momentum,
the continuous rotating telescopic-pair phase and signed-work refinement,
frame/link-order/endpoint covariance, hidden-spin refusal under a large common
boost, a real 32-row driven graph, polynomial
collapse/representation/resource refusals and every cancellation checkpoint.
This is an ideal bilateral length source: finite motor/storage/efficiency laws,
unilateral cables, contact/impact coupling, exact energy preservation,
symplecticity and validated whole-path state enclosures remain outside this API.

### Spherical and mixed closures between independent native bodies

`rigid_rattle::distance::graph::advance_joint_graph` accepts `GraphJoint::Spherical`
and `GraphJoint::Distance` in one coupled graph. A `SphericalJoint` makes two
body-fixed material points coincide while leaving relative attitude free.
It contributes three Cartesian closure rows, never a zero-length radial rod
with an undefined normal. Each body keeps its finite mass and full inertia;
shared-body responses couple all three axes and every attached rod through the
same scaled mobility, fs-la factorization, local correction/backtracking and
native impulse/free-step/impulse engine used by the static and driven entries.
Existing distance-only and driven APIs retain their original numerical path.

`JointGraphStep` contains the original native graph receipt and an authored-joint
to scalar-row map. Both reaction inventories preserve actual equal-and-opposite
body-point impulse events in x/y/z order for each ball joint, including every
shared-body work cross term. No positions, velocities or energy are repaired.
The mixed receipt separately pairs signed constraint work from double-double
relative COM and spin components before their point velocities are rounded.
It also returns the total requested angular impulse: exact coincident endpoints
have no internal couple, but finite admitted closure residuals can leave a small
couple that must stay visible. Original native actual momentum changes and work
residuals remain available; a small closure residual is not an exact angular
momentum or energy guarantee.

Initial and final spherical gap/speed admission uses the whole Euclidean vector
norm, as does the per-event spherical impulse cap. These physical bounds are
independent of the selected world axes. Cartesian gaps retain
`(COM_second-COM_first)+arm_second-arm_first` and rates retain relative COM
velocity plus both native spin terms, so a distant origin or common boost cannot
hide local geometry or velocity before admission. The state still uses native
binary64 coordinates and momenta, with their original representability limits.

The existing maximum of 32 independent scalar rows admits, for example, ten
ball joints plus two rods. Rank, body count, trial, physical branch and byte
limits remain authoritative. One fixed-size stack row bank and one fallibly
reserved joint map extend the original reused buffers; there is no new numerical
engine or per-trial allocation. The reported logical payload includes the map.
Every late refusal or cancellation publishes no body or partial joint account.
This bounded entry supports smooth regular bilateral spherical/static-distance
graphs. Driven spherical anchors, mixed actuator commands, welds, redundant
rows, contact/impact/friction coupling, exact energy preservation, symplecticity
and global or validated whole-path certificates remain outside it.

`tests/rigid_joint_graph.rs` checks an independent three-axis unequal-mass
ball-chain load-sharing solution without welding spins; a rotating finite-body
ball-pair continuous solution with second-order refinement and reversal;
200 full-tensor mixed-graph steps with native P/L/work accounts; frame, boost,
body/joint order and endpoint covariance; exact distance-only/free-path parity;
local geometry and hidden-spin admission at a `2^54` origin or common boost;
whole-vector caps, actual connected 32-row admission and rank/resource refusals;
and cancellation at every observed successful transaction checkpoint.

## Articulated point/plane impact response

`articulated_impact::resolve_articulated_plane_impact` resolves one declared
link-fixed point against one prescribed plane at a supplied configuration. It
normalizes the finite nonzero outward world normal, computes point velocity
and the full coupled normal effective mass, and consumes the single-contact
Newton law owned by `fs-contact::impact`. The inverse-mass operation reuses the
existing zero-bias free-base ABA path; no dense generalized mass matrix is
assembled. A unit-response probe is not a physical impact and is not compared
with endpoint joint-speed or actuator-force limits. The final physical impulse
still passes the existing atomic impulse owner's endpoint admission.

The borrowed input state and all positions remain unchanged. A strictly
negative gap refuses; a gap above the declared activation tolerance or a
nonapproaching contact receives zero impulse. Returned state includes base and
joint recoil, actual endpoint relative normal speed, world linear/angular
impulse, and before/after kinetic energy. The prescribed plane supplies work
`W_plane = p*n.dot(v_plane)`; the receipt exposes `Delta K - W_plane + D`, where
`D` is the contact law's nonnegative dissipation. This accounts for an actively
moving plane without misclassifying its work as spurious energy creation.

This is one frictionless impact event, not collision-time detection, position
repair, a persistent-contact force, a simultaneous-impact solver, or a
Moreau–Jean time stepper. Six independent impact regressions exercise offset
COM/effective inertia, rotated frames, articulated recoil, elastic/inelastic
and intermediate restitution, moving-plane work, activation/penetration, and
actual endpoint limits. Existing articulated impulse and multibody unit tests
also cover the shared inverse-mass extraction. The broad contact and multibody
beads remain open for their remaining dynamics and constraint deliverables.

## Simultaneous articulated contact response

`articulated_contacts::resolve_articulated_contacts` advances the velocities of
up to 128 free articulated systems through one event with at most 64 declared
frictionless normal contacts. A contact joins a link-fixed point to a prescribed
plane or another dynamic link point, including two links in the same tree.
The adapter normalizes each world normal and checks dynamic witnesses share a
normal line of action so paired impulses cannot invent an angular couple.
Strict negative gaps refuse; activation tolerance admits proximity only.

The contact Delassus matrix is assembled through the shared zero-bias ABA
inverse-mass owner, probing only the one or two systems involved in each
contact. A single globally admitted contact-law solve produces all impulses.
The adapter accumulates equal/opposite link loads before applying one physical
impulse update per system; it never performs sequential pairwise state updates
or assembles a dense generalized-coordinate mass matrix. Existing joint-speed
admission applies to the actual endpoints, not the unit-response probes.

Before returning any states, the adapter measures endpoint contact velocities
from the actual updated articulated kinematics and checks them against the
coupled predictions. Global linear and angular momentum changes must equal the
prescribed-plane impulses. Plane work is signed and explicit; measured total
kinetic-energy change minus that work plus the contact loss must close within
the declared numerical energy allowance and roundoff. A failed body update,
contact-law refusal, unresolved balance, or cancellation publishes no candidate
state; input positions and borrowed velocities remain unchanged.

The roundoff allowance scales with actual event storage increments, moving-plane
and prescribed-support work, constitutive storage change and physical loss. It
does not scale with unchanged absolute kinetic inventories. A massive coasting
spectator cannot excuse an incorrect contact loss. Original full-inertia energy
and support work remain retained until the assembled comparison; a small resolved
transfer inside a tree survives cancellation with another tree's larger signed
transfer. This changes event admission, without inferring heat from a numerical
work defect. Focused `articulated_contacts::momentum_tests` cases exercise both
of these boundaries through the native response and simultaneous solve.

Independent analytic regressions cover two and three coupled contacts on one
body, an unequal-mass three-body impact chain, an internal prismatic-link pair,
contact-order invariance, moving-support work, world-frame/common-velocity
transforms, and a rebound candidate that satisfies local constraints but would
inject global energy. The latter refuses while its explicitly requested
inelastic counterpart succeeds. Redundant directions remain a typed numerical
refusal. No friction, collision-time search, position repair, continuous contact
force or nonsmooth timestep claim is made.

### Mixed free and prescribed-base normal impacts

The additive `articulated_contacts::mixed::resolve_mixed_articulated_contacts`
entry admits simultaneous normal contacts between independently free trees,
stationary supported mechanisms and externally driven mechanisms. Each
`MixedArticulatedContactSystem` declares its actual free or prescribed base;
the existing native ABA inverse-mass maps supply its contact mobility. All
contacts share the original geometry, simultaneous normal law, accumulated
link impulses and endpoint admission, including its redundant unique-velocity
response with representative impulses and admitted total work. The original
free-only entry retains its native response type through the result's default type
parameter.

Prescribed bases retain their actual native reaction force and couple, in
body and world coordinates, through `MixedArticulatedImpulseResponse`.
Their support work is reported separately from plane work and mechanical
impact loss. The shared global momentum residual subtracts both plane and
base impulses; its energy residual is actual kinetic change minus both work
sources plus admitted loss. A moving support can supply energy without
being mistaken for passive-contact energy creation. No fictitious support
mass, independent link mass, sequential collision response or post-event
velocity projection is introduced.

The original limits of 64 contacts and 128 systems, full-tensor and
offset-COM treatment, actual endpoint speed limits, nonnegative-gap admission,
global restitution passivity and cancellation remain active. Locked or
insufficiently resolved active normal mobilities refuse through the original
solver. This entry adds instantaneous normal momentum exchange; the adjacent
circular and directional mixed-boundary entries extend it to friction.
Collision detection and sustained-contact trajectories remain separate owners.

Focused native tests cover free-only response parity, a free body striking a
supported hinge with physical support recoil, two coupled contacts sharing
one nondiagonal mobility, a moving and rotating full-tensor support under
frame/origin/boost changes, and atomic rank, speed-limit and cancellation
refusals.

### Shared material-point contact geometry

`articulated::contact_geometry` owns link-fixed witnesses, positive separating
world normals, point impulse directions, signed gap/support speed and the
common-line check for dynamic witness pairs. It consumes original admitted
`Kinematics`; `GeneralizedKinematics::native_kinematics` exposes that same
snapshot to bounded generalized-inertia consumers. Native
`articulated_contacts` reexports the existing public contact types and calls
this owner for its normal preparation and shared point/normal helpers. Thus
the dry ABA adapter and a fluid-augmented generalized adapter use the same
material points and sign conventions without duplicating contact geometry.
Strict negative gaps and invalid or noncollinear witnesses retain their native
refusals. The existing seven normal-contact regressions continue to cover
coupled body/self contacts, moving support, frame transformations, global
energy admission and atomic refusal through this shared geometry path.

`prepare_contact_frame` supplies the same normal/tangent/normal-cross-tangent
rows to the native Coulomb adapter and the fluid-augmented consumer. It
normalizes and orthogonalizes the original tangent hint, retains both actual
material witnesses, and measures all three incoming relative speeds from
native link twists. Borrowed kinematics slices admit either owned maps or
references without cloning the trees. The original
`ArticulatedFrictionalContact` type and native public entry points remain in
place. Finite-gap dynamic tangential couples are not hidden by witness
projection; the consuming mechanics must retain its world-angular-momentum
audit.

### Coupled articulated Coulomb contact

`articulated_contacts::friction::resolve_articulated_frictional_contacts`
extends the same native mechanics to simultaneous normal and two-tangent
impulses. It consumes the original `ArticulatedContactSystem`,
`LinkContactPoint` and `ContactCounterpart` types. Stable distinct contact IDs
canonicalize the full contact matrix, solver input and accumulated link loads.
The normal-only API and friction API share their actual link Jacobians,
zero-bias ABA inverse-mass actions, physical impulse application and independent
endpoint energy/momentum audit. Joint/base recoil, rotated inertia and all
normal/tangent and cross-contact couplings remain in `J M^-1 J^T`; no lumped
link mass or sequential friction clipping replaces them.

The adapter calls the existing `fs-contact::frictional_multicontact` owner
without changing its non-associated Coulomb law, restitution/target semantics,
cone and complementarity gates, or determinacy evidence. Returned velocities
are measured again from the actual articulated endpoints. All systems are
updated privately and published together only after physical momentum, signed
support work and total kinetic loss close. Unit probes do not consume joint
speed limits; the actual endpoint retains that admission. Dynamic witnesses
are never moved to create a common point, so an activated finite-gap friction
couple can still fail the angular-momentum gate.

There are at most 21 contacts and 128 systems. Each contact row probes only its
one or two involved trees with linear-time ABA; the bounded contact matrix is
dense, but the generalized articulated mass matrix is never assembled.
Velocity-only determinacy retains individual impulses as an admissible
representative. Its support work requires a common translation among active
prescribed supports; independently moving supports refuse without an additional
work-determinacy witness. No local heat or pressure authority is inferred.

G1/G3/G4 regressions in `tests/articulated_contacts/friction.rs` cover the
zero-friction normal-owner limit, analytically coupled off-centre stick/slide,
prismatic and revolute recoil, a shared-body conveyor under contact permutation
and tangent-frame rotation, two dynamic counterparts, and late cancellation
with deterministic retry. This is a frozen-configuration point impulse with one
common static/kinetic coefficient. It adds no contact discovery, time stepper,
velocity-dependent friction, finite-patch spin loss, or heat allocation.

### Stateful articulated shear and reversible contact storage

`articulated_contacts::tangential::advance_articulated_tangential` consumes the
shared `fs_tribo::TangentialLaw` and its retained checkpoint on one fixed physical
tangent. The original articulated ABA supplies the exact projected mobility
`W = J M^-1 J^T`; the frozen-configuration midpoint force compliance is
`dt^2 W / 2`. The constitutive owner's work-conjugate average traction applies
`dt F` to that same mechanics. Joint/base recoil and both dynamic counterparts
remain in the response; no lumped link mass or scalar approximation to a 2-D
mobility is introduced. A translating prescribed support supplies signed work.

Actual endpoint velocities reconstruct the midpoint slip after native impulse
application. Contact storage before/after, constitutive work, physical heat,
support work and native kinetic-energy change close one account. The shared
mechanical audit now includes explicit reversible contact-storage change, so a
preloaded contact can accelerate resting bodies by releasing stored energy.
Normal-only and ideal Coulomb impact calls retain zero contact-storage change
and their existing passivity semantics. Lost mechanical energy from an impact
is not automatically labeled friction heat by the new consumer.

Input bodies, law and history are borrowed. One complete trial returns both the
constitutive checkpoint and all native body responses; failure and cancellation
publish neither. Geometry, force, time, slip-speed, compatibility and work budgets
are explicit. Foreign law/source histories retain their original refusals;
transverse force or motion refuses because only one coordinate is assembled.
The caller keeps the physical tangent/history association unchanged and publishes
all returned state together. This operation provides no pose integration,
normal-response constraint, changing-load history transfer, multiple coupled
shear ports, finite-patch spin/rolling, thermal allocation or contact discovery.
Dry sticking is step-averaged; time refinement must resolve branch transitions
and intrastep reversals independently of work closure.

Six G1/G3/G4 groups in `tests/articulated_contacts/tangential.rs` cover an
independently derived off-centre elastic load line, preload energy release,
prismatic/revolute joint recoil, dynamic-pair momentum and frame covariance,
the analytic Mindlin backbone and retained reversal memory, dry sliding heat,
and late refusal/cancellation with exact retry. All six and the thirteen prior
normal/Coulomb groups pass natively using the original mechanics, contact and
numerical solver bodies. This focused source-selected validation excludes
unrelated runtime consumers and does not establish a full Cargo/DSR gate.

### Continuous articulated dry friction and exact arrest

`articulated_contacts::tangential::flow::advance_articulated_dry_friction`
resolves up to six simultaneous scalar guide/contact ports with the original
`fs-contact` continuous static/kinetic owner. Source-declared Coulomb laws and
normal loads supply both capacities; complete native `J M^-1 J^T` retains joint,
rotational and shared-body recoil. Exact zero-crossing events end sliding
segments. Resting coordinates resolve their coupled static reactions together,
including breakaway under held drives, with the original ambiguity refusal.
There is no capture-speed threshold or midpoint velocity reflection at arrest.

Continuous signed travel, unsigned distance, per-port friction work and held
actuator work survive every event. Integrated drive and friction impulses apply
once to the original ABA mechanics. Its existing momentum/energy audit now
accepts explicit source work without relabeling it as negative contact storage
or negative heat; other event consumers retain exactly zero source work.
Prescribed-plane work remains a separate physical boundary account. Canonical
contact IDs make caller permutation deterministic, and cancellation, dependent
mobility or a physical/work limit returns no partial native candidate.

The source law, normal loads, body/joint geometry and projection are held for
the phase. Returned travel is not an SE(3) pose update; changing geometry and
loads still require outer refinement. The scalar ports do not approximate a
two-axis Coulomb disk, and no normal constraint or thermal partition is inferred.
Analytic regressions in `tests/articulated_contacts/tangential/flow.rs` exercise
off-centre arrest, driven rest versus tiny slip, coupled static recoil,
finite-body and moving-support work, original prismatic/revolute response,
deterministic retry and atomic refusal.

### Two-axis stateful articulated shear

`articulated_contacts::tangential::planar::advance_articulated_planar_tangential`
extends retained shear to both axes of the existing validated `TangentFrame`.
The original ABA assembles the full two-by-two `J M^-1 J^T`, including both
cross terms, and the shared `fs-tribo` matrix embedding resolves the original
constitutive law against `dt^2 W / 2`. Both impulses are accumulated before one
native mechanical update. Dynamic counterpart and same-tree joint recoil stay
in that response; a transverse stored force is neither dropped nor separately
integrated. No generalized mass matrix is assembled.

The original endpoint velocities independently reconstruct both actual slip
components. One audit closes contact storage, physical heat, total articulated
kinetic energy, prescribed-support work and world linear/angular momentum.
The bounded matrix solve polls cancellation at constitutive evaluations, and
the caller receives body responses and retained history only as one complete
candidate. Original material-law refusals and structured matrix-solve failures
remain visible. The caller must keep the physical frame/history association;
this operation does not invent frame transport, normal evolution, a pose step,
multiple-contact coupling or finite-patch rolling/spin. Scalar Mindlin retains
its transverse-motion refusal. Dry sticking remains a midpoint-slip condition.

Eight G1/G3/G4 groups in `tests/articulated_contacts/planar.rs` cover independent
off-centre two-axis load lines, elastic preload release, nonproportional layer
plasticity, moving support work, dynamic-pair covariance, an independent
Coulomb multiplier oracle, rank-one internal-joint mobility, frictionless
motion and unchanged retry after late refusal/cancellation.
All eight and the nineteen prior contact groups pass in the focused native
original-source build; this does not establish a full Cargo/DSR gate.

### Simultaneous retained shear at multiple articulated contacts

`articulated_contacts::tangential::planar::multiple::advance_articulated_tangential_contacts`
advances an explicitly bounded list of original contact/law/history pairs in
one common physical interval. The native ABA assembles every entry of the full
`2N`-row `J M^-1 J^T`; cross-contact blocks remain in `C = dt^2 W / 2`.
The shared prepared tribology owner resolves all retained constitutive histories
against that operator. Only then are all impulses accumulated into one native
update per articulated system. Contact count and `(2N)^2` storage are admitted
before probing or allocating the contact matrices; no generalized mass matrix
or independent-contact approximation is introduced.

Each actual endpoint supplies its own midpoint slip and storage/work/heat check.
The total physical heat and storage change enter the original global kinetic,
support-work and world momentum audit exactly once. Every candidate history and
all body responses return together, in the caller's original order. One source,
geometry, history, solve, actual-work or cancellation failure publishes none.
The original per-site normal load, temperature and pressure remain explicit;
different durations refuse. Normal response, pose and physical frame transport,
interface ownership, discovery, rolling/spin and thermal allocation remain with
their existing owners. Dry sticking remains a midpoint-slip condition.

Six G1/G3/G4 groups in `tests/articulated_contacts/planar/multiple.rs` cover
opposing and three orthogonal shared-body interfaces against an independent
dense elastic solution, contact permutation, dynamic-counterpart momentum,
retained nonproportional plastic histories, manufactured simultaneous stick/slide
and unchanged whole-system retry after late refusal or cancellation.
All six and the twenty-seven prior contact groups pass in the focused native
original-source build; this does not establish a full Cargo/DSR workspace gate.

### Objective transport of retained articulated shear

`articulated_contacts::tangential::planar::transport` adds framed single and
simultaneous contact proposals. Each site borrows an accepted
`FramedTangentialState`, declares its endpoint contact frame and supplies an
integrated **common material spin** angle. The shared original `fs-tribo`
transport rotates retained elastic, plastic and dry accumulated-slip vectors
isometrically, preserving native state ownership and scalar history. A tangent
gauge change does not rotate the physical preload; common spin is not relative
spin and adds no slip, traction, work or heat by itself. Antipodal normal changes
refuse instead of selecting an arbitrary material rotation.

After transport, the unchanged original single or simultaneous articulated
proposal assembles its full ABA contact mobility and applies one accumulated
impulse. The result moves each original candidate history into its endpoint
frame once. All bodies and framed histories return together; refusal or
cancellation at any stage leaves the accepted inputs unchanged. Count/matrix
caps precede transport allocation. Per-contact and global actual-state work
checks include transport energy roundoff under the existing mechanical energy
allowance, even when the separate transport allowance is larger.

This operation is an explicit transport-then-shear split at supplied frozen
endpoint geometry. The shared minimum-normal rotation followed by common spin
is an endpoint path approximation; endpoint normals and a scalar spin integral
do not reconstruct arbitrary noncommuting frame motion. The enclosing owner
supplies material/contact identity, pose/normal evolution, the spin integral
and refinement across changing frames or constitutive reversals. No contact
discovery, normal solve, relative-spin friction, rolling law or thermal store
is introduced, and scalar Mindlin gains no unsupported two-axis transport.

Six G1/G3/G4 groups in `tests/articulated_contacts/planar/transport.rs` exercise
world-force invariance under tangent gauge changes with plastic memory,
independent tilted-preload recoil, common-spin and three-dimensional world
covariance, zero fictitious spin slip/heat, distinct simultaneous frame
transports against a dense elastic oracle, and exact retry after frame/source
or late mechanical/cancellation refusal.
All six and the thirty-three prior contact groups pass in the focused native
original-source build; this does not establish a full Cargo/DSR workspace gate.

### Reusable articulated response for coupled endpoint trials

`PreparedArticulatedTangentialContacts` retains the original immutable native
contexts, admitted contact rows, full ABA mobility and shared constitutive
factor for one frozen interval. Repeated `advance` calls repeat neither the
inverse-mass probes nor matrix admission/factorization. They still invoke the
original accumulated native impulse application and all actual endpoint,
storage/heat/support-work and world momentum checks on private pristine
scratch. A failed or cancelled trial cannot contaminate the next proposal.
Frames, support motion, body pose/velocity and duration belong to preparation;
changing them requires a new owner. Model and joint arrays keep their original
immutable borrows, while native base contexts are copied snapshots.

`PreparedArticulatedFramedTangentialContacts` additionally transports each
original accepted history once before repeated environment trials. No trial
ages history or reapplies common material spin. The original native law must
accept that checkpoint with unchanged stored energy; parameter work and source
continuation require their own explicit owner. Original transport residuals
remain in every candidate's accepted-state energy audit. Prepared results
borrow both dense contact matrices; `into_owned` copies them only when the
caller retains a final report. The existing one-shot functions keep their
owned-vector API and move the prepared matrices into their results.

Three G1/G3/G4 groups in `tests/articulated_contacts/planar/prepared.rs` cover
independent full elastic recoil, unchanged repeated histories and matrices,
exactly-once native frame transport, snapshot and final-report ownership, and
exact retry after late speed, temperature, shape/interval or cancellation
refusal. This preparation alone neither solves temperature nor advances time.
All three groups and the forty-four prior groups pass together in the focused
native original-source target; this is not a full Cargo/DSR workspace gate.

## Native rigid-body finite-patch tangential impulses

`tangential_friction::apply_rigid_friction` binds supplied positive normal
patches to original COM/body-frame `RigidBodyState` and `MassProperties`
values. It derives `B = J M^(-1/2)`, including translational mass,
orientation-dependent rotational inertia, shear lever arms, torsional couples,
and all shared-body coupling. No authored effective mass substitutes for body
inertia. Relative point velocity, relative spin, and common spin come from native
point kinematics; the existing unregularized tangential island determines the
impulses.

Actual native point-impulse updates and equal/opposite torsional couples are
staged on copied bodies. Independent native endpoint-velocity, kinetic-energy,
midpoint-work, total linear/angular momentum, and whole friction-energy checks
gate one atomic body/history publication. Separated contact points refuse unless
within the explicit geometry tolerance; any resulting angular-momentum defect
still independently refuses. Changing-coefficient work remains an external signed
port, while only physical microslip loss is friction heat.

Poses stay fixed: this is an impulse phase requiring the surrounding integrator's
geometry advance and normal law, not a complete contact time integrator,
restitution solver, or Moreau–Jean stepper. The broad nonsmooth normal/impact
deliverables remain open. Three physical regression groups cover rotated
anisotropic off-center bodies with analytic shear/torsional recoil, multiple
contacts sharing a body with exact input-order invariance, and late balance,
cancellation, duplicate-work, and geometry refusal without body/history mutation.

## Event-resolved sphere/plane trajectories

`sphere_plane::advance_sphere_plane` advances the original momentum-form rigid body through ballistic detection, restitution, frictional impact, persistent normal support, driven rolling/sliding and liftoff. It admits one isotropic sphere with physically enclosed COM inertia, a fixed-normal plane moving at constant velocity, a constant world force at the centre, and a common static/kinetic ideal Coulomb coefficient. Timing and impulse response remain with the original `fs-contact`/`fs-tribo` owners.

`advance_sphere_plane_with_support_friction` admits the original `fs-tribo::FrictionLaw::Coulomb` with distinct static and kinetic coefficients for sustained contact. Its validated range is `0 <= mu_k <= mu_s`; source calibration remains the caller's responsibility. Zero-slip rolling uses the static limit, while already moving contact uses the kinetic limit until resolved arrest. Collinear arrest then re-establishes rolling inside the static limit or starts reverse sliding above it. This includes zero kinetic friction with nonzero static capture. The scene's original coefficient still controls instantaneous impacts; a support law does not silently redefine impact response. Velocity-dependent laws require another evolution owner and refuse here.

During support, only the normal applied force is balanced. With zero applied couple and tangent inverse mass `kappa = 1/m + radius^2/I`, rolling requires traction `-(F_t/m)/kappa` inside the circular friction limit. Collinear sliding retains analytic constant-force motion, exact slip-arrest times and resolved reversal. Frictionless support admits arbitrary tangent drive.

Resolved changing-direction sliding uses the original isotropic Coulomb disk in a velocity-implicit discretization, bounded by `max_sliding_step_s` and `max_sliding_substeps`. Actual translation, world momentum, contact-point slip and torque include the original spin inertia. Each accepted step uses its resolved constant force in the native trajectory; endpoint Coulomb feasibility does not assert an exact continuous force direction. The curved trajectory is first order and arrest is resolved at a timestep boundary. Refinement is necessary. Coarse trials that would make friction assist the initial slip are halved, including apparent endpoint arrest. Supercritical transverse drive cannot be silently converted to persistent rolling. Rejected trials count toward the work cap.

For distinct coefficients, the transverse noncapture condition retains the kinetic limit: driving at or above that limit cannot acquire a finite-time zero-slip endpoint merely because static friction could hold a separately supplied resting state. Thus the hysteresis band admits both persistent rolling and continuing slip according to actual history. The exact critical transverse case also retains its physical noncapture condition. The five focused groups in `tests/sphere_plane/support_friction.rs` check analytic breakaway, arrest/re-stick/reversal and zero-kinetic capture; unchanged impact response; curved noncapture; and moving-support work, split continuation and atomic refusals.

All five and the fifteen prior sphere/plane groups pass with the pinned compiler against the original selected native dynamics/contact sources and full real tribology/math dependencies. This is focused native verification; it does not establish the full Cargo/DSR workspace gate.

Both routes evaluate lost work from actual midpoint relative slip on each constant-force segment, never from the fictitious impulse-law loss at the force-driven free velocity. Curved-step work belongs to the declared time discretization and converges with its sliding timestep. Rolling traction has zero relative loss; the moving plane supplies signed work. No mechanical loss is automatically assigned to heat.

Attitude uses an exact exponential for commuting spin/torque and bounded second-order Lie midpoint substeps otherwise. Newly supplied penetration refuses. Each curved supported segment checks its complete quadratic normal-gap path against the original floating geometry allowance, in addition to actual endpoint kinematics. Private continuation retains only bounded contact roundoff; no pose or momentum is projected.

Every accepted advance independently checks total linear impulse, orbital-plus-spin angular impulse about the world origin, and `K_after - K_before = external_force_work + external_torque_work + plane_work - dissipated_energy`. Cancellation, unresolved geometry/work, nonadvancing clocks and exhausted event, attitude or sliding budgets publish no partial state. Accumulating rebounds do not alter restitution.

Fifteen focused native integration groups pass: the eleven original ballistic, support, collinear, incline, spin, frame, attitude and refusal groups plus four curved-flow groups. An independent separable continuum solution for `u' = a - B*u/|u|` verifies velocity, position, spin and loss convergence at 16/32/64 steps, with error ratios 1.8–2.2; a separate 32/64/128 ladder verifies arrest followed by driven rolling. Coarse-step passivity, supercritical non-arrest, moving-plane work, split continuation, cancellation and late sliding-budget rollback are also exercised. These are focused actual-source native tests, not a full Cargo/DSR workspace gate. General multicontact Moreau-Jean dynamics, time-varying or body-held torque, anisotropic bodies, multiple planes, finite-mass plane recoil, rolling/torsional resistance, arbitrary-shape certified CCD, source calibration and contact adjoints remain outside this owner.

### Applied world torque during sphere/plane contact

`sphere_plane::advance_sphere_plane_with_torque` and
`advance_sphere_plane_with_torque_and_support_friction` add a prescribed
constant world couple about the sphere COM to the original event engine.
Its angular acceleration contributes `(torque/I) cross arm` to the original
contact-point driving acceleration. The existing Coulomb owner therefore
resolves motor-driven static rolling, traction-limited breakaway, braking,
slip arrest and reversal, including distinct static/kinetic support limits.
Torque-updated pre-impact spin enters the original impulse law. Translation,
geometry, restitution and contact reaction ownership remain unchanged.

Each accepted segment integrates isotropic world spin under the sum of
applied and contact torque. Signed couple work uses the shared compensated
midpoint pairing of actual segment spin rates and couple impulse. The public
receipt separates this work and direct angular impulse from centre-force
work, plane work and contact dissipation. The complete orbital-plus-spin and
kinetic/work admissions include all these channels. Zero applied couple
retains the original arithmetic. Commuting attitude is exact; changing-axis
attitude and curved sliding retain their separate declared refinement bounds.
Cancellation and budget refusals return no partial state.

Five focused physical groups cover powered rolling on a moving conveyor and
spatial covariance, normal-axis spin without fictitious point-contact heat,
traction breakaway/static hysteresis/reversal, torque-dependent impact friction
and split continuation, transverse sliding against independent Newton/Euler RK4,
and zero-couple parity with late refusal. The ideal spherical point contact
does not acquire torsional friction, finite-patch rolling resistance or a
time-varying/body-held motor law from this extension. Material-derived rolling
retains its existing zero-applied-couple API.

### Motor torque through directional support and impact

`sphere_plane::directional::advance_sphere_plane_directional_with_torque`
carries the same prescribed COM world couple through original directional
support and optional instantaneous-impact callbacks. The native
`directional::impact::advance_sphere_plane_directional_impacts_with_torque`
supplies the existing distinct sustained and impulse laws. Actual driven point
acceleration and pre-impact spin determine traction and impact response.
The fixed plane texture, source-phase/native-slide/attitude bounds and atomic
history publication retain their original owners. Returning no custom impulse
preserves the scene's circular impact law. Existing entries delegate zero couple.

Force work, signed couple work, prescribed-plane work, sustained slip heat and
impact loss remain distinct. Four G1/G3/G4 groups in
`tests/sphere_plane/directional/torque.rs` cover weak/strong texture breakaway,
zero point-slip drive with nonzero force/couple work, torque-dependent impact
with frame/conveyor covariance, and late source/cancellation/forged-work/
nonfinite-torque refusal. Native rustfmt parses the changed sources and
independent analytic ledgers close; native fs-mbd test execution remains unrun
under the unavailable current dependency closure. No Cargo/DSR gate is claimed.
The constant world-couple and isotropic point-contact model boundaries remain.

## Native rigid-body simultaneous frictional impacts

`frictional_impact::apply_rigid_frictional_impacts` binds the simultaneous
normal/Coulomb owner to original rigid-body states. Actual masses, full body
inertias, orientations and COM-relative contact arms determine the full
normal/tangent Delassus matrix. Assembly contracts only shared-body influence
blocks, avoiding products over unrelated body coordinates. Stable contact IDs
also canonicalize native impulse accumulation.

`PointKinematics::relative_velocity_world` and
`velocity_relative_to_world` reconstruct relative motion from the original
COM velocities, angular velocities and point arms. Double-double products and
reductions cancel common translation before rounding the final relative
velocity. Both native impact admission and endpoint verification use this
query, including directional impacts and material-source callbacks. Absolute
point-velocity diagnostics can round away a small rotational contribution;
they are not used to decide whether that contact is approaching or slipping.
If the resulting required momentum jump cannot be represented, the original
native impulse and endpoint gates still refuse the complete event.

Dynamic pairs receive equal/opposite original point impulses; translating
prescribed planes supply an explicit external impulse, angular impulse and
signed support work. Independently reconstructed native endpoint velocities,
kinetic/midpoint work, total linear/angular momentum and whole-event energy
must agree before any body is published. Initial penetration, unresolved
contact or any late balance failure leaves all bodies unchanged. Distinct
dynamic witnesses are not silently moved to a fabricated common point; a
tangential couple due to an admitted finite gap can still fail angular balance.

The shared circular/directional impact assembly reduces original fixed-inertia
kinetic products across every body before rounding. Moving-support work and
relative-contact midpoint work retain independent product layers through the
complete energy gates; scalar per-body impulse reports do not rebuild storage.
`RigidKineticEnergyChange::finish_retained` exposes the same retained reduction
to coupled finite-material/heat owners. Representable small physical loss
survives opposing large transfers, while actual momentum rounding and source
work errors remain rejection defects. Focused regression definitions cover a
20.25 J pair loss that separately rounded body reports change to 20.5 J, and a
10.125 J moving-support loss hidden by scalar storage/work subtraction. The
small production arithmetic tests do not execute the complete impact solver.

Four native physical tests cover off-center rotated anisotropic inertia with
independently derived normal/tangent cross-coupling, sticking/restitution and
moving-support work, a simultaneous three-body chain with exact input-order
invariance, and penetration/redundancy/cancellation/late angular refusal.
The contact core's SPD and determined-response scope remains in force.
This point-impact phase holds poses fixed, requires supplied geometry and
coefficients, and makes no automatic finite-patch preload, torsional-friction,
collision-discovery or impact-to-heat allocation claim.

### Independent native rigid contact islands

`frictional_impact::islands::apply_rigid_contact_islands` partitions supplied
contacts by shared dynamic bodies and applies the original simultaneous
normal/Coulomb solver to each connected component. Its directional counterpart,
`apply_rigid_elliptic_contact_islands`, uses the same decomposition and the
existing friction ellipse owner. A prescribed support has no reciprocal state,
so a common floor or conveyor does not connect independent dynamic bodies.
Contacts sharing any dynamic body remain one coupled solve, including supplied
inactive witnesses. An oversized connected component refuses before solving;
it is never split into sequential pair collisions.

The scene body/contact caps bound graph storage and the number of solves;
`contact.maximum_contacts` is the dense owner's per-island limit, defaulting
to 21. An explicitly larger contact cap, matching scalar matrix budget and
Krylov work budget can admit up to 256 contacts through the original coupled
Newton/Coulomb law. More than 63 active scalar rows require positive-definite
mobility and the original weighted contraction witness; bounded small singular
mobility and velocity-only paths are not silently enlarged. Graph storage is
linear in scene bodies and contacts. Dense
mobility storage belongs only to individual islands, not a square matrix for
the whole scene. Stable contact IDs determine island and global output order,
with original body/contact index mappings and every native determinacy receipt
retained. Redundant impulse representatives gain no local-load or heat authority.

All native body states remain private until every island and the original
absolute whole-scene momentum/work gates pass. The numerical allowances do not
grow with island count. Compensated reductions preserve the summed physical
increments, source/support work, contact losses and residuals. A late solver,
balance, index, capacity or cancellation failure publishes no partial update.
`step_rigid_contact_islands` and
`step::elliptic::step_rigid_elliptic_contact_islands` reuse the original force
kick, native drift, mandatory whole-path callback and complete-step balance
gates. This does not introduce a new contact law or relax any connected
solver's applicability.

### Directional simultaneous rigid impacts

`frictional_impact::elliptic::apply_rigid_elliptic_impacts` consumes a principal
coefficient pair and tangent axis at every original contact. It shares the
original body Gram assembly, canonical impulse application and complete
endpoint/momentum/work admission with the circular entry point. Its receipt
keeps physical impulses and velocities separate from the contact solver's
explicitly transformed determinacy witness. An invertible tangent congruence
preserves unique-impulse versus unique-velocity authority; redundant individual
loads still do not become unique pressure or heat allocations. The transformed
contact owner's explicit work-range witness also establishes unique total
relative work. Compatible rotating supports, independent islands and exact
frictionless tangent motion remain admissible; incompatible redundant work
returns the original transformed `UnresolvedWork` refusal. Individual support
loads and source-work shares remain representatives.

`apply_rigid_elliptic_impacts_with` resolves each geometrically active source
once, in canonical contact-ID order, from the original preimpact point velocity,
exact assembled normal/tangent frame and both actual witness positions. The
query retains the original input index. Inactive candidates are unqueried and
retain their histories, although their endpoint velocities still include the
other contacts' impulses. Source candidate states and original body momenta
publish together only after the entire event passes; cancellation, source
failure, solver limits and late angular/work refusal discard all candidates.
No event duration, normal pressure, finite patch or frictional heat fraction is
inferred. Coefficients are explicitly sampled before impact, not solved as an
inner velocity-dependent impact law. A single zero principal coefficient is
an unsupported singular ellipse; two zero coefficients are frictionless.

Four focused native G1/G3/G4 groups cover a rotated ellipse with full body
inertia and off-centre normal/tangent coupling against an independent endpoint
law; moving-support/frame covariance; a simultaneous three-body chain and
circular-limit identity; redundant directional support/work authority; and
source-index/kinematics, inactive-history, cancellation/budget and late-audit
rollback followed by retry. These remain frozen-configuration point events,
without collision discovery, persistent directional time stepping or automatic
impact-loss disposition.

## Configuration-admitted native rigid contact steps

`frictional_impact::step::step_rigid_frictional_contacts` advances actual
native momenta and poses through a first-order constant-wrench/gravity kick,
simultaneous admitted normal/Coulomb impulses, and the original zero-force
energy/world-momentum midpoint drift. The step independently measures initial
and endpoint signed gaps, advances translating support planes, and refuses
every negative supplied gap without projection.

A mandatory geometry-owner callback must admit the complete proposed native
path and feature coverage. Its native partial-drift sampler is a tool, not a
CCD certificate; initial/endpoint checks or finitely many samples alone do not
establish continuous nonintersection. Callback refusal, unresolved contact,
drift/work failure and cancellation leave every body unchanged.

Receipts separate actual kick work, constant-force displacement work, their
signed first-order splitting difference, moving-support work, unallocated
impact loss, drift defects, and independent whole-step energy/momentum balances.
Neither event loss nor integration loss is relabeled as physical friction heat.

Four native trajectory tests cover a finite-height block with a supplied
admissible center of pressure at incline stick/threshold/slip conditions,
sustained nonpenetrating support, three-rung first-order trajectory/work/loss
refinement, original asymmetric free drift, moving supports and late atomic
refusal. The independent sliding fixture gives `x(1)=1+dt`, kick work
`6+44dt`, and event loss `4+44dt` for `dt=1/8,1/16,1/32`.
The geometry callback bounds the entire cuboid path rather than checking only
the supplied point.

This is a caller-admitted first-order composition within the contact owner's
SPD/determinacy domain. Automatic manifold/center-of-pressure discovery,
arbitrary CCD/Painleve resolution and a complete Moreau–Jean integrator remain
separate work. Geometry and time refinement remain explicit requirements.

### Variable-mass rocket translation

`rocket::RocketState::advance` integrates a point vehicle with explicit dry
and propellant masses, engine mass flow, world effective-exhaust vector,
gravity acceleration and other external force. The effective exhaust vector
includes nozzle pressure thrust once. Piecewise-constant forcing uses the
Tsiolkovsky logarithm and its separate analytic position integral; external
force is integrated against the changing mass, rather than divided by an
average mass. A bounded polynomial series preserves tiny burned fractions.
Propellant exhaustion splits a step into powered motion and coast at dry
mass, with zero engine thrust afterward. The operation is deterministic and
allocation-free, validates representability and polls cancellation before
publication. G1/G3 tests cover independent analytic velocity/position,
constant-force integration, burnout, split steps, tiny burns, frame changes
and Galilean invariance. This point-vehicle operation supplies no atmosphere or aerodynamic
coefficients and no accuracy guarantee for forcing that changes within a step.
The separate `rocket::rigid` extension below adds native attitude and staging
under explicit mass-distribution and exhaust-momentum assumptions.

## Native finite-patch Hertz impulse mechanics

`frictional_impact::hertz_patch::apply_rigid_hertz_impulse` applies the original
`fs-tribo` Hertz tangent/spin impulse law to one finite rigid body and a second
finite body or prescribed moving patch. The supplied compressive normal impulse
sets the constitutive budget and is not applied again. The caller supplies the
fixed radius, pressure-impulse shape, coefficient and coincident material
witnesses. Geometry, source motion and every accepted native endpoint retain
their original dimensional gates.

The complete radius-scaled three-coordinate mobility is assembled from the
original rigid point-response factors and central spin couple, including full
inertia tensors and off-centre translation/spin coupling. Each touched body
receives one native resultant impulse containing both the material-point moment
and free spin couple. The original full normal cross response and actual
before/after normal velocities remain visible; this conditional phase does not
impose a restitution target or solve coupled normal contact.

`PreparedRigidHertzImpulse` caches the fixed native snapshot and mobility for
repeated normal-budget/coefficient trials. Changing the radius requires new
preparation. Only the one or two touched bodies are staged, and successful
preparation/solve/publication uses fixed stack storage. An allocation is used
only when retaining a failed audit report.

Signed prescribed-support work includes both patch-centre translation and
normal-axis rotation. Actual native momentum increments determine kinetic change
without subtracting large endpoint energy totals. Native relative midpoint
work, constitutive loss, support work, endpoint response and linear/angular
momentum balances are checked independently. Cancellation or any failed gate
leaves all caller bodies and reusable preparation unchanged. Mechanical loss
remains unallocated; contact duration, elastic-memory evolution, heat partition
and subsequent pose integration belong to their existing owners.

## Coupled normal and finite-patch impacts on native bodies

`frictional_impact::hertz_patch::impact::apply_rigid_hertz_impact` connects the
coupled normal/Hertz owner to the original finite body pair or prescribed
moving/spinning patch. `PreparedRigidHertzImpulse::solve_impact` reuses a
fixed native snapshot for coefficient/restitution trials. The cached full
four-coordinate Gram response contains normal, both tangent directions and
radius-scaled spin, including all original inertia and off-centre cross terms.

The normal impulse is determined together with sliding/spin and applied
exactly once. Each touched body receives one summed native point impulse and
central spin couple. The conditional and coupled APIs share the original
native application, polarized momentum/kinetic account and support-work
calculation. Normal restitution/complementarity, actual tangent/spin
endpoints, complete linear/angular momentum, total relative work and the
constitutive energy result are independently gated before publication. Moving
support work includes its normal and tangent translation plus angular motion.

The native report exposes the applied normal impulse. Its conditional patch
receipt keeps its own conditional work, while the native and coupled totals
refer to the complete event. The coupled elastic limit retains signed
roundoff within the explicit energy allowance; genuine active impact candidates
refuse. No negative loss is hidden by clipping or relabeled thermal energy.
Existing conditional tangent/spin calls still apply no normal impulse and
retain their original strict passivity gate.

Only one or two native states are staged, poses stay fixed, and every failure
or cancellation preserves the caller's bodies. Prepared geometry uses fixed
storage; the coupled contact solve additionally owns its bounded original
numerical factors. Both layers retain the same patch impulse bounds. This
single-impact model assumes a caller-supplied effective Hertz radius and
pressure-impulse shape; a transient footprint, compliant impact duration,
multicontact restitution and heat allocation remain separate physical owners.

## Complete finite-mass two-sphere trajectories

`sphere_pair::advance_sphere_pair` advances two original rigid-body states
through their earliest sphere contact and the entire remaining duration. Both
bodies have finite, possibly unequal masses/radii, physically enclosed isotropic
inertia and arbitrary spin. Their applied acceleration is the same constant
world vector. Relative flight is linear, so the shared swept-sphere detector
and a nonclosing accepted impact require at most one event and two analytic
flight phases. Torque-free attitudes use exact exponentials without substeps.

The adapter supplies both bodies' translational and rotational effective mass
to the existing normal and Coulomb owners. Equal and opposite impulses act at
the actual facing surface points. Actual endpoint velocities, relative impulse
work, total linear momentum, orbital-plus-spin angular momentum, applied-force
work, mechanical loss and final separation are independently admitted before a
complete result is returned. New initial overlap refuses; positions are never
snapped together, and cancellation returns no partial state.

Focused physical tests cover unequal-mass restitution, two-tangent spin/friction
recoil, interior crossings, noncollision, common-acceleration and frame/velocity
covariance, and atomic refusals. This is a floating analytic two-body lane with
equal static/kinetic Coulomb coefficients, not a many-body contact network,
unequal-acceleration detector, compliant/adhesive model, anisotropic solver,
rolling-resistance model or automatic heat allocation.


### Applied torque through two-sphere flight and impact

`sphere_pair::advance_sphere_pair_with_torques` adds a prescribed constant
world couple about each sphere's centre of mass. Isotropic spin angular
momentum advances analytically on each free-flight segment. The resulting
pre-impact spin enters the original normal/Coulomb response, and the
post-impact spin starts the remaining torque interval. Signed applied-torque
work (including braking) and the direct couple impulse are reported separately
from centre-force work and orbital moments. Complete kinetic/work and
orbital-plus-spin angular-momentum balances include both flight segments and
the original impact loss; no state or residual is projected.

Commuting angular velocities use the exact integrated rotation. Changing spin
axes use native quaternion exponentials at angular-velocity midpoints, with
second-order attitude accuracy and explicit `SpherePairTorqueConfig` controls.
The total noncommuting substep budget covers both bodies and both flight
segments. Conservation does not certify phase accuracy: callers must resolve
attitude with the declared time bound. Budget exhaustion, invalid input and
late cancellation return no partial state. The original zero-torque entry
retains its existing flight arithmetic.

Focused regressions cover analytic spin-up, braking and reversal; torque-driven
pre-impact slip and frictional recoil with simultaneous force work; frame
covariance; changing-axis attitude convergence against an independent RK4
reference; and budget/cancellation rollback with deterministic retry. This
extension retains isotropic enclosed inertia and common constant acceleration.
State-dependent or body-held couples, anisotropic inertia, sustained contact
and many-body torque-driven trajectories remain outside this entry.

### Native redundant support and prescribed-support work

`frictional_impact::apply_rigid_frictional_impacts` consumes the contact core's explicit impulse/velocity determinacy receipt. Its original shared-body Gram mobility makes a unique contact velocity image sufficient for unique native body momenta, even when individual contact impulses are indeterminate. Such impulses remain an admissible representative; the report does not promote them to unique local pressure or heat.

A velocity-only event carries the contact core's explicit impulse-work range witness, proving that its total relative work is independent of the selected contact-impulse representative. The original Gram mobility then also determines total prescribed-support work from the actual kinetic change. This admits compatible translating or rotating supports, mixed finite-body contacts and independent islands without requiring a common support translation. Exactly frictionless tangent coordinates can be excluded from the work proof, so distinct harmless tangent belt speeds remain admissible. Incompatible redundant source work returns `Contact(UnresolvedWork)`; the directional owner retains the original transformed error. Individual support loads and source-work allocations remain representatives. Existing actual endpoint-velocity, kinetic/work, linear/angular-momentum and cancellation gates still precede atomic body publication.

Native G1/G3 groups in `tests/frictional_impact/redundant.rs` cover a finite-height, three-point body in stick and gross slide with independently derived load transfer, energy and momenta; frictionless duplicate normal rows with restitution; conveyor work; and a rotating platform coupled to a finite-body collision under frame changes and contact permutations. Directional regressions retain the transformed witness and circular-limit behavior. Generic rank/feedback, configuration, contact discovery and physical heat-partition boundaries remain explicit.


### Automatic native sphere collision trajectories

`frictional_impact::spheres::advance_sphere_contacts` advances complete COM-centred native spheres under explicit uniform `Gravity`. `advance_forced_sphere_contacts` additionally admits one finite constant world force at each original centre of mass, in body-slice order. Every pair has an immutable coefficient/restitution law; stable sphere IDs canonicalize discovery and simultaneous impulse order. Equal force/mass accelerations retain the original straight relative-path query. Different accelerations use `fs-contact::swept_sphere::accelerated::sphere_path` to isolate every interior distance extremum and bracket first entry on each complete quadratic relative path.

No pose or witness is projected. Up to 64 outward floating-time neighbours resolve nonnegative actual witness gaps AND complete path gaps; unresolved geometry refuses. Original radial material arms retain mass, principal inertia and orientation. Nearby contacts enter one simultaneous normal/Coulomb event at the declared spatial resolution. Native free rotation keeps its original approximation and acceptance gates.

Gravity and additional centre-force work use actual accepted displacements. Linear and world angular impulses remain separate. Gravity torque integrates the original total mass moment, including total applied centre force; each body-specific force torque integrates its actual piecewise quadratic path across collision-induced velocity changes. Each segment and the full trajectory close kinetic work, linear momentum and orbital-plus-spin angular momentum before atomic publication. Impact loss remains unallocated mechanical dissipation.

Four new G1/G3/G4 groups extend ten original collision/gravity cases: independently forced elastic entry with analytic gravity/force work and angular impulses; manufactured curved 3-D frictional entry with an independent original-inertia Coulomb impulse and post-impact trajectory; accelerated departure/return, permutation and split continuation; and late event-budget/cancellation rollback, deterministic retry and refusal of unresolved loaded resting pairs. These focused original-source native tests do not establish a Cargo/DSR workspace gate.

This bounded floating analytic lane retains all-pairs discovery. It does not claim arbitrary shapes, interval-certified CCD, force fields, time-varying loads, body-held couples, persistent supports, rolling resistance, thermal allocation or full Moreau-Jean dynamics. Zero-speed contact requiring a sustained reaction refuses rather than integrating into overlap. The broader `frankensim-ext-contact-nonsmooth-lane-oh0i` remains open.

### Applied world couples in automatic sphere collision trajectories

`frictional_impact::spheres::advance_wrenched_sphere_contacts` accepts one
held world couple about each sphere COM alongside the existing body-specific
centre forces and common gravity. Prepared native `PreparedSpatialLoads`
and `RigidBodyIntegrator::step_spatial` own the full-inertia endpoint torque
kicks, gyroscopic rotation and couple work. Rotational substeps advance a pure
spin factor; the original native constant-acceleration centre drift runs once
per flight segment. Refining torque therefore does not perturb the analytic
centre path used by the existing accelerated CCD queries. Torque-updated spin
reaches the original simultaneous normal/Coulomb event, and post-impact spin
starts the remaining forced flight.

`SphereTorqueControls` supplies separate time and whole-call work bounds.
The budget counts all bodies, event segments and outward root-bracket trials;
zero couples retain original free-drift arithmetic and consume no torque
steps. Angular phase and asymmetric-body energy accuracy still require
refinement. Signed free-couple work and direct angular impulse remain separate
from gravity/centre-force work, orbital force moments and unallocated impact
loss. Native drift defect stays numerical. Both segment and full-trajectory
work/momentum checks precede publication; late budget or cancellation refusal
leaves every input body unchanged.

Four physical groups cover analytic simultaneous three-sphere motor/friction
coupling with gravity and centre-force work, body permutation and world-frame
covariance, independent full-inertia Euler/quaternion refinement, and braking,
zero-couple parity and atomic refusal. The API holds world couples over the
declared interval; body-held or state-dependent torque, sustained supports,
rolling resistance and thermal allocation remain outside this event lane.


### Sustained redundant support with the exact contact linearization

The native rigid consumer exercises the analytic semismooth contact action and cached numerical range preconditioner without changing force/contact/drift physics or acceptance tolerances. A three-foot body now resolves sustained normal-load transfer and friction below, exactly at, and above the Coulomb threshold without a supplied centre of pressure or numerical impulse drift. The complete containing-prism path, original native mass/inertia, actual momenta and work/energy balances remain checked.

`g1_sustained_tripod_support_resolves_load_transfer_without_a_supplied_pressure_centre` advances three load regimes at 8, 16 and 32 steps over one second: 168 accepted steps. For the sliding regime its independent first-order discrete reference is x(1)=a(1+dt)/2. Position errors against constant-acceleration continuum motion are 0.00625, 0.003125 and 0.0015625 m, halving on each refinement; subthreshold and threshold runs retain sticking within the declared gates. The containing prism admits the original COM-centred mass distribution and inertia, and the contact points remain on its actual lower face.

All eight native impact tests, six automatic sphere-trajectory tests and four original contact-step tests pass against the final core, with no filtered cases or relaxed gates. These are source-selected native checks with the repository-pinned nightly; they do not establish a full Cargo/DSR workspace gate, a general Moreau-Jean scheme or unrestricted contact determinacy.

### Native four-foot support with indeterminate pressure

The original native impact and force/contact/drift owners now consume the contact core's global candidate-direction admission for determined body motion with redundant frictionful normal reactions. A uniform 2 m cube with four corner supports retains every contact and its actual mass/inertia. No normal load distribution or centre of pressure is supplied, and individual foot reactions remain an admissible representative.

Three new native groups verify sticking, exact Coulomb threshold and sliding; oblique two-axis circular friction and both pressure moments; canonical contact permutation; and sustained load/trajectory refinement. The latter advances loads below, at and above threshold at 8, 16 and 32 steps over one second (168 accepted steps), checking whole-cube path separation, actual body momentum and work. Sliding position errors are 0.00625, 0.003125 and 0.0015625 m against the continuum reference, halving with the original first-order time step.

All eleven native impact groups, all ten gravity-aware sphere groups and all four existing contact-step groups pass against the frozen new core. The old duplicate-contact refusal fixture now uses different friction coefficients on the unresolved shared reaction, preserving a genuine ambiguous-motion refusal instead of rejecting determined motion merely because impulses are redundant. These focused source-selected native checks do not establish a full Cargo/DSR gate. Existing support-work restrictions, fixed contact-law scope and unallocated impact-loss semantics remain in force.


## Native rigid rocket attitude and conservative staging

`rocket::rigid` composes the original analytic variable-mass translation with
native `UnitQuaternion`, `MassProperties`, `RigidBodyState` and the existing
momentum-preserving free-rotation map. Its mass model has a fixed common COM,
shared principal axes, and inertia linear in remaining propellant. Dry and
propellant principal moments are checked for physical admissibility.

Body-directed equivalent exhaust velocity supplies total thrust, including
pressure thrust once. Its application arm supplies the matching thrust torque.
A separately declared diagonal co-rotating exhaust moment and optional swirl
account for exported spin; neither is silently inferred from the thrust arm
or reservoir inertia. Symmetric torque/outflow kicks surround free rotation
at midpoint inertia, and the original analytic translation uses midpoint
attitude. The smooth coupled map is second order, bounded by caller-owned
substeps, allocation-free, and segmented at actual burnout. Results retain
scalar axial impulse separately from directed world impulse, expelled-mass
closure and world angular-momentum accounting. Cancellation publishes nothing.

`rocket::rigid::staging` returns both separated rigid bodies. It admits the
specified child masses, COMs and full parallel-axis inertia against the parent,
applies equal-and-opposite internal recoil, and reports actual mass, linear
momentum, orbital-plus-spin angular momentum and separation-work residuals.
Children retain the declared common material body axes; their principal axes
may differ when using the full-tensor model. Incompatible geometry or
nonrepresentable bookkeeping is refused.

### Full-tensor rocket stage separation

The existing `RocketStage<I>` and `StageSeparation<I>` now retain the same
sealed `RocketInertiaModel` as their parent. The diagonal model remains the
default; `TensorRocketInertia` stages preserve independently oriented dry and
propellant tensors without projecting them onto common principal axes.
Admission compares every entry of the complete parallel-axis tensor, as well
as separate dry/propellant inventories and the common COM. Both returned
children keep their original inertia laws and continue through the native
full-tensor depletion and flight owner.

The original assembly release still owns the inherited rigid velocity field,
common-point equal/opposite recoil, spin plus orbital momentum, and separation
work. The adapter remains allocation-free for two stages. Three focused
G1/G3/G4 regression groups cover independent tensor/impulse/work calculations,
no-recoil release followed by powered flight of both descendants, an admissible
off-diagonal-only inertia mismatch, and late cancellation with deterministic
retry. These are numerical fixtures, not calibration or a resolved explosive
separation model; geometry/contact admission remains the caller's operation.

These operations do not infer internal relative-flow momentum/Coriolis
transport, slosh, shifting COM/principal axes or arbitrary outlet inertia
tensors. Prescribed equivalent thrust is not a complete internal-fluid model,
and constant inputs over one step do not imply exact varying-force flight.

Native verification with the pinned compiler used the actual dynamics core
and actual `fs_math` dependency. All 17 focused rocket tests passed (nine new
rigid/staging tests and eight prior translation tests). Checks include analytic
translation/burnout, free asymmetric rotation, co-rotating mass-loss spin,
thrust-arm torque, measured second-order position/velocity/attitude refinement,
frame and Galilean invariance, staging momentum/work and atomic refusals.

## Automatic forced sphere scenes and translating supports

`frictional_impact::sphere_scene::advance_sphere_scene` discovers contacts and
advances complete COM-centred native spheres against other spheres and
prescribed infinite planes with constant translation. It composes the shared
bounded spatial broad phase with the original native force/contact/drift step
using independent dynamic contact islands. Each body retains its actual
mass, orientation, spin and physically enclosed principal inertia. Applied
world force, body torque and uniform gravity remain explicit inputs.

Stable sphere and plane IDs identify the material-law lookup. A missing law
refuses the update; it cannot disable a discovered collision. Contact arms are
the actual radial sphere surface arms. The normal restitution target uses
pre-kick incoming velocity, so a resting body under gravity requests support
without converting its force kick into repeated restitution bounces. The
original simultaneous Coulomb solver still admits actual endpoint velocities,
passivity and support-work determinacy.

Every proposed drift checks the complete straight relative sphere path and
the affine sphere/plane gap, including the actual native endpoints. A rejected
path can suggest an earlier encounter time; every shorter retry repeats the
force kick, contact solve and complete-path admission because the post-kick
velocity depends on the new duration. Initial overlap refuses. No position
projection repairs penetration. Nearby contacts enter the same simultaneous
solve within the declared activation resolution.

The complete requested duration is staged. Cancellation, missing laws,
geometry failure, exhausted comparison/candidate/plane-query/step/retry
budgets, or failed native conservation gates leave caller bodies unchanged.
The report retains accepted native step receipts, geometry work counts,
contact identities and whole-trajectory linear/angular momentum and energy
accounts. Force impulse work, force-displacement work, splitting defect,
torque work, prescribed-support work, contact loss and native drift defect
remain separate. Contact loss is unallocated mechanical dissipation, including
possible normal timestep loss; it is not automatically friction heat.

Whole-trajectory conservation uses compensated sums of the original native
kick/contact/drift increment residuals, with one unchanged global tolerance.
Absolute kinetic energies and orbital momenta do not define this balance by
subtraction: a large common translation or a distant origin cannot erase a
small contact loss or support torque. Source/support work, support impulses
and loss diagnostics also use compensated chronological reductions. Multiple
locally admitted impulse errors can still exhaust the whole-call allowance
and refuse the complete trajectory without changing caller-owned bodies.

The caller's maximum step controls the original first-order force/contact/drift
approximation. Geometry retries do not establish an adaptive time-error
estimate. This numerical lane supports straight drift of spheres and fixed
normal translating planes; it does not claim interval-certified CCD, arbitrary
shapes, rotating supports, deformability, adhesion, rolling/torsional resistance
or thermal allocation. The original simultaneous solver's coupled-contact cap
applies to each connected island. Explicit contact, scalar matrix and Krylov
budgets enable larger positive-definite connected groups; the default remains
21. The separately declared scene contact cap also bounds independent contacts.
Accepted native step reports retain
the complete per-island solver and body receipts. The broader contact and
multibody beads remain open.

Six focused native G1/G3/G4 groups cover resting support with nonzero material
restitution, falling contact and settling, two-layer stack load transfer,
independently derived incline rolling/sliding trajectories, conveyor impulse
and work, hidden interior collisions with encounter-time retries, and complete
rollback after retry/plane-query exhaustion or cancellation. All six pass with
the original native dynamics, contact, Newton/Krylov, geometry and math owners;
the seven shared broad-phase groups pass separately. This selected-source
native verification uses the pinned nightly. It does not establish a full
Cargo/DSR workspace gate; the required sibling constellation was unavailable.

### Directional friction in automatic sphere scenes

`advance_elliptic_sphere_scene` uses the same discovery, retry, path admission
and whole-call staging as `advance_sphere_scene`, dispatching the original
simultaneous elliptic contact-island solver. `SphereSceneEllipticLaw` supplies
two principal Coulomb coefficients and an explicit
`SphereSceneTextureDirection`. The texture can be held in the world frame or
attached to either contacting sphere by stable source ID and body direction.
The actual staged attitude transports attached texture at each accepted
interval; rejected retries cannot advance it. Tangent projection uses the
actual contact normal. Zero, nonfinite, normal-parallel or unrelated texture
directions refuse instead of inventing a material axis.

Both coefficients define one coupled maximum-dissipation ellipse. A zero
principal coefficient retains exact free sliding in that direction; equal
coefficients recover the circular law. This does not introduce separate
static/kinetic phase history. Normal restitution still uses pre-kick incoming
velocity. Every contact retains native finite-body recoil, moving-support work,
complete island coupling and the original momentum/energy admission.

The material frame is held during each native first-order interval and refreshed
after accepted drift. The same authored time resolution and straight-path
geometry limits apply. Mechanical contact loss remains unallocated.
Generic report/error defaults preserve the scalar API; elliptic aliases expose
the original directional receipts. `tests/sphere_scene/elliptic.rs` covers
circular parity, a frictionless principal direction, rotated moving supports,
attitude transport, sphere-pair recoil and complete rollback.

## Sustained loaded articulated contact

`articulated_contacts::step::step_articulated_frictional_contacts` adds a
first-order force/contact/drift composition on the original free-floating
articulated states. Held internal joint efforts, link-frame wrenches, and
uniform gravity produce one coupled ABA impulse response at the actual
configuration. Gravity acts at each original link COM, including its moment
about the link origin. The original force/state validator checks joint effort
limits before conversion to impulses. The simultaneous normal/circular-Coulomb
owner then resolves every supplied contact together; the native unforced Lie
midpoint owner advances every base and joint with those accepted velocities.
Internal actuator and dynamic-contact reactions retain native whole-tree recoil.

The actual initial and endpoint signed gaps must be nonnegative. A mandatory
geometry callback admits the entire candidate drift path, feature coverage,
and continued applicability of the frozen material witnesses and normals.
It may query the original numerical path; finite sampling and endpoint tests
alone are not CCD. Prescribed planes translate once at their declared velocity.
No position or velocity projection repairs overlap, late path refusal, contact
failure, or cancellation. All models and input states are borrowed; publication
of every endpoint and the caller's geometry receipt is atomic.

Reports retain the original load/contact responses, exact kick work, signed
support work, native unallocated contact loss, and independently measured final
linear/orbital-plus-spin angular balances. Exact held-joint displacement work
and uniform-gravity COM-displacement work expose their separate first-order
splitting defects. The gravity angular splitting term compares the initial kick
with its exact integral along the accepted linear whole-tree COM drift. A
rotating body-frame wrench's continuous path work is not inferred from endpoint
poses. Native contact loss includes the force/contact splitting loss when a
support removes the free load kick; it is not automatically continuum friction
heat. Unforced drift error remains separately measured. Independent user gates
bound per-tree/aggregate drift energy and complete-step work/momentum error.

The lane admits at most 128 trees, 4096 total links, and the native 21-contact
cap, with bounded cancellation polls and O(links) tree passes. It adds sustained
resting/sliding contact for caller-admitted paths, not automatic arbitrary-shape
collision discovery, a general Moreau--Jean integrator, exact continuous-force
motion, or heat partition. Numerical activation resolution remains explicit,
and time-step refinement is necessary. Both broad contact and MBD beads remain
open. Focused G1/G3/G4 tests in `tests/articulated_contact_step.rs` exercise a
supported actuated carriage, sliding with independently conserved carriage
momentum and first-order trajectory refinement, multi-tree load transfer,
moving-support work, contact permutation, full-scene rotation, native effort
limits, and late geometry/cancellation/penetration refusal.

## Event-aware articulated trajectories

`articulated_contacts::trajectory::advance_articulated_events` composes the
original free-floating Lie-midpoint motion owner with the original simultaneous
Coulomb contact owner. All free bases and scalar joints advance under declared
held loads to each supplied event time, receive one simultaneous impulse solve
at fixed poses, and continue with the accepted native velocities. Equal-time
sets must be combined into one event; sequential resolution of simultaneous
contacts is refused. Event-free intervals are uniformly subdivided under an
explicit maximum smooth duration and interval budget. An event at the start or
end is admitted. A zero-event, single-interval call reproduces the original
smooth endpoint exactly.

The returned states can begin the next window. Aggregate physical momentum is
recomputed from actual link states about a caller-chosen fixed world origin;
angular accounting includes spin and the moving COM's orbital contribution.
Native external impulse, actuator work, gravity/body-wrench work, prescribed
support work, total impact loss and smooth numerical work defects remain
separate. Dynamic contact pairs cancel globally. Native event receipts retain
contact determinacy and passivity evidence without inventing a heat partition.
The shared contact admission uses uncancelled original link linear, spin and
orbital contributions at both endpoints for its existing roundoff multiplier:
a reaction wheel's opposite base/rotor momenta must not collapse that scale to
zero and reject a physically valid centered impact.

The window admits at most 128 systems, 4096 total links, 64 events and 4096
common smooth intervals; each event retains the native 21-contact bound. No
generalized dense mass matrix is assembled. Caller state is immutable, and all
stages and native event results remain private until the final cancellation
poll succeeds. A late invalid event, physical admission failure or cancellation
returns no partially advanced state. G1/G3/G4 tests in
`tests/articulated_trajectory.rs` cover analytic gravity rebound, oblique
friction/spin, two-tree recoil, an actuated reaction wheel, moving-support work,
multiple-event continuation/replay, and late refusal/cancellation.

This owner consumes caller-located times and geometric witnesses; it does not
establish no missed collisions, exact event times, persistent contact, a
nonpenetrating between-event path, or full constrained-integrator accuracy.
Moving-plane positions are supplied at each event and are not advanced twice.
Smooth step size still requires refinement. Native world-origin angular
impulses are translated to the selected origin, so distant-origin floating
cancellation remains a numerical limitation rather than a hidden projection.

## Automatic spherical articulated carriers against translating planes

`articulated_contacts::trajectory::sphere_planes::advance_carrier_planes`
advances one declared spherical carrier with its original internal articulated
rotors through automatically discovered plane impacts. Every link origin and
centre of mass must coincide with the base origin; only fixed/revolute joints
are admitted. Link mass second moments must fit the declared spherical radius.
That consistency check does not prove enclosure of an unknown shape: the caller
declares the sphere as the exposed surface enclosing all internal links. External
link forces are refused; internal actuator efforts and body torque loads remain
supported. Uniform gravity then supplies the whole-system COM acceleration.

Each original Lie-midpoint step is limited by the earliest closing root from
the native ballistic gap owner over all supplied planes. Every accepted segment
checks the complete quadratic sphere/plane gap, including interior minima, and
the actual native endpoint. Fixed-normal planes translate once from their
declared window-start position. Isolated grazing needs no impulse; an initially
departing contact can return. Floating roots whose actual surface witnesses
overlap retry up to 64 earlier floating times, with no position/velocity repair.
The observation is floating analytic geometry, not an interval CCD certificate.

At an event the actual radial root-link material arms are rebuilt from native
orientation. Nearby planes enter one simultaneous original circular Coulomb
solve in stable-ID order. That solve retains shell spin, joint/base recoil,
original inertia, global passivity and signed translating-support work. All
smooth actuator/external work, impact loss, numerical work defects and spatial
momentum accounts remain separate in the original trajectory receipt. Spatial
activation resolution permits early impulses only at nonnegative gaps within
that resolution; no penetrating configuration becomes admissible by tolerance.

The complete requested window is private until the final cancellation poll.
At most 21 planes, 4096 links, 4096 attempted native steps and 64 impacts are
admitted, with caller-owned geometry and event budgets. Failed/retried native
steps count toward the work cap. Loaded slow contact, unresolved same-time
events, exhausted event budgets (including chatter/Zeno), impossible geometry,
and cancellation return no partially advanced state. Sustained support requires
a different contact lane. There is no claim for other-body collision coverage,
protruding internals, rotating planes, rolling constraints, certified rotational
error, or general Moreau-Jean dynamics. Both broad contact and multibody beads
remain open.

Focused G1/G3/G4 cases in `tests/articulated_carrier_planes.rs` cover analytic
gravity rebound, a hidden interior collision between separated endpoints,
motor-driven rotor recoil with frictional shell spin, simultaneous corner
contact and plane permutation, translating-conveyor work, repeated rebound and
split continuation, full-scene rotation with a body torque, isolated
grazing/departure/return, and atomic budget, cancellation and geometry/model
refusals. Native selected-source checks are
separate from the unavailable full Cargo/DSR workspace gate.

## Material-derived two-axis rolling with finite-body recoil

`tangential_friction::rolling::apply_rigid_rolling` consumes the original
`RollingLossLaw::ViscoelasticSphere`, normal patch, ordered interface, static
`FrictionLaw`, work key, native rolling checkpoint, rigid-body states and mass
properties. Both bodies recoil. Their actual arms, orientations and principal
inertias determine the full four-coordinate shear/rolling mobility. Eliminating
the two no-slip reactions gives the full two-axis rolling Schur response;
projected Gram factors avoid cancellation from subtracting inertial matrices.
The isotropic sphere law acts along the actual tangent relative angular
velocity. Neither rolling speed nor an effective inertia is supplied by callers.

Original Hertz/material admission checks load, area, radius, temperature,
frequency and the quasistatic small parameters. The actual sphere arm must equal
its radius along the contact normal, and both contact points must coincide
within the explicit geometry gate. A loaded contact requires zero initial and
final point slip/separation and relative normal-axis spin within the dimensional
gates. The native static capacity bounds the required traction throughout each
held interval, including interior extrema. An infeasible reaction refuses instead of forcing
rolling or substituting an independently saturated friction component.

The material coefficient remains owned by the original card. A bounded scaled
two-dimensional rotation resolves the native symmetric constrained response;
the two orthogonal rates each decay exactly with their original eigenvalue.
No eigenvalue-gap division or near-isotropic cutoff changes the physical model.
The original material owner integrates both components in one receipt: signed
rate integrals supply the couple impulse, while squared-rate integrals supply
physical heat and one cumulative-state/checkpoint update. The full vector speed
must satisfy the original source applicability gates.

Anisotropic reaction components can cancel initially and peak later. Their
squared-norm derivative reduces to a quadratic in the relative exponential,
so both endpoints and every interior stationary traction are checked. Native
normal/spin mobility also determines the continuous omitted-coordinate paths;
their possible interior extremum is checked before publication. These are
floating-point analytic evaluations, not interval-certified bounds.
The user retains the time-resolution ceiling and explicit substep work cap,
but damping stiffness introduces no extra steps. Viscous arrest remains
asymptotic; zero rate, zero viscosity and zero load give exact zero loss without
snapping or a velocity regularizer. Only the sphere generates bulk heat. The
report exports that heat and its original body index without creating a thermal
reservoir or partitioning it into the rigid plane.

Actual native point impulses and equal/opposite couples independently close
kinetic change, midpoint work, physical heat, zero static-reaction work, and
total linear/orbital-plus-spin angular momentum. Candidate bodies and native
history publish together after the final cancellation poll. Late balance,
source, mode, traction and budget refusal leaves both unchanged.
Only the two contacted native bodies are staged and audited in fixed-size
storage. Unrelated scene bodies are neither scanned nor copied, and their large
kinetic energies cannot erase the pair's physical loss in an energy subtraction.
Report kinetic and momentum totals belong to that pair; their changes equal
the whole-system changes. Original body indices and source roles are preserved.

The rolling, joint Hertz/rolling, resolved collinear and shared-island reports
retain `kinetic_energy_change_j` from the actual accepted momentum endpoints.
The native fixed-inertia polarized pairing accumulates all participating bodies'
compensated products before rounding, so unchanged common motion and opposite
large transfers cannot erase smaller resolved kinetic changes. Absolute totals
remain descriptive outputs; work/heat gates and the original finite-caloric
rolling adapters consume the retained increment. Commanded work and material
heat remain independent, and rounded-away momentum still enters the defect.

This is one held-geometry, held-normal-load pure-rolling phase. The surrounding
owner still advances poses and normal contact. It does not solve contact
discovery, lift-off, externally driven no-slip constraints, simultaneous
microslip/spin, multiple patches, or arbitrary-speed viscoelastic rolling. A
second shear reaction cannot be superposed on this same phase. Frozen-geometry
decay is exact for isotropic and anisotropic response; this does not certify a
complete trajectory. Native regressions compare both body momenta and material
heat to independent constrained mechanics above the midpoint sign-reversal
threshold with a one-step work budget, and across split caller intervals.
G1/G3/G4 groups in `tests/rolling_friction.rs` use independently eliminated
principal-direction mechanics for both-body recoil, three interval partitions,
nearly repeated eigenvalues, zero components, nonreversing decay, exact zero limits, world and
surface-order covariance, source/traction/mode refusal, and loaded checkpoint,
late audit, cancellation and exact retry.
The tests also require isolated-pair recoil, heat and history despite an energetic
untouched body and nondefault contact indices, and reject a physical two-body
case whose interior traction exceeds capacity although both endpoints fit.
These focused original-source checks are separate from a full Cargo/DSR gate.

## Material-derived supported sphere rolling

`sphere_plane::rolling` advances the original supported sphere pose and momenta
under the existing material-derived Kelvin–Voigt rolling moment. The original
Hertz patch, source density/radius, held temperature and frequency, and
quasistatic validity remain binding. The constant-load vector equation
`(m + I/R²) du/dt = F_t + (T × n)/R - (A N/R) u` is integrated analytically,
including driven turning and reversal. The original entry point retains zero
applied torque; `advance_viscoelastic_sphere_rolling_with_torque` admits a held
world motor or braking torque. Both endpoint static-traction norms bound the
entire exponential path; infeasible rolling refuses rather than projecting
away slip. Native linear and world-origin angular impulses, actual force work,
signed motor work, moving-plane work and continuous bulk torque work are audited.
The nonnegative heat integral uses exponential moments, not an energy defect.
Normal torque accelerates normal spin and its work is retained separately;
normal spin contributes no deformation heat in this rolling law.
Changing-axis attitude alone uses bounded Lie midpoint steps;
commuting attitude is exact. This is a homogeneous-sphere, held-load no-slip
phase with asymptotic coast-down, not sliding, impact, finite-time viscous
arrest, torsional loss or thermal transport. Generated heat belongs to the
card's deforming sphere and is already included in the native advance loss.
Focused physical groups cover exact coast-down, motor startup, braking and
reversal, driven work, normal-spin work without spurious heat, world/boost
covariance, continuation, elastic rigid-body limits, source/traction refusal,
late cancellation and attitude refinement, including motor startup from rest
whose rolling and normal-spin accelerations evolve differently.

`advance_viscoelastic_sphere_rolling_until_breakaway` returns the complete
requested interval or the audited physical prefix at its first static-traction
exit. The traction path `m*a0*exp(-rate*t)-F_t` is a line segment in force space;
the original scalar ballistic-contact owner locates its disk-boundary root.
Capacity-scaled line distance preserves early exits, and a reverse query from
terminal traction preserves late exits without subtracting nearly unit decay
fractions. A boundary state that initially unloads keeps its later outward
crossing; terminal traction on the disk is approached asymptotically. Event
time is independent of trajectory timesteps. Shared `expm1` evaluation retains
small early traction changes under nearly cancelling motor and force drives.
The accepted endpoint stays on the feasible side of floating-point root
evaluation; this is not an interval-certified enclosure.

The result retains all accepted native work, material heat, source identities,
and the remaining horizon. Breakaway adds a bounded, zero-impulse native
`Sliding` event and supplies required static traction, capacity and instantaneous
material rolling torque for handoff. Immediate breakaway preserves pose and
momenta exactly, with zero duration, impulse, work and heat. Source validity
still bounds the whole prefix, and event/attitude budget exhaustion or late
cancellation publishes nothing. The original full-interval APIs retain strict
traction refusal. No post-breakaway sliding is propagated: its deformation and
surface-friction losses need a separately admitted owner, and the material
torque is never silently discarded. Focused cases cover exact force- and
motor-driven prefixes, immediate and inward-boundary events, early/late clocks,
force scaling, continuation, moving-frame work, source limits and cancellation.

## Native rigid-body directional friction impulses

`tangential_friction::directional::{apply_rigid_directional,
apply_rigid_directional_with}` applies a held-load, frozen-geometry
directional friction phase to two original finite rigid bodies. Their
masses, full inertia tensors, poses and center-of-mass-relative contact
arms determine the full tangential mobility `W = J M^-1 J^T`; the midpoint
constitutive compliance is `C = dt^2 W / 2`. No scalar averaging of the
off-diagonal response is performed.

The native entry point retains the original elliptical law and history.
The callback entry point admits a retained source owner using the same
history, released motion, matrix and selected branch. A source that owns
speed/path coordinates must close them against actual solved motion
inside the callback. Exactly zero initial tangential speed selects the
static-or-sliding branch; moving contact selects kinetic continuation.

The actual opposing impulses update both bodies. Independently
reconstructed endpoint kinematics, per-body and total momentum balances,
midpoint work and contact storage/dissipation gate publication. Kinetic
change is computed from momentum increments and mean velocities rather
than subtraction of large absolute energies. The report exposes normal
recoil for the surrounding normal owner. Exactly two bodies and accepted
history are staged and published together after cancellation and all
physical gates; failure leaves them unchanged.

A phase that assists endpoint slip or produces moving midpoint stiction
refuses as an unresolved slip event. Exact resolved stopping is admitted.
The trajectory owner must refine or locate arrest; this adapter performs
no event search, pose integration, normal impulse solve, changing-load
integration, spin-friction solve or thermal-store update.

`tests/rigid_directional.rs` checks full-tensor off-center response against
independent unit-impulse columns, equivalent body-frame parameterizations,
work/heat balance, exact stopping under large common translation, overstep
refusal and atomic callback/cancellation/geometry failures.

## Generic retained tangential laws on finite rigid bodies

`tangential_friction::apply_rigid_tangential` couples an original `TangentialLaw`
and its accepted checkpoint to the complete native two-axis point mobility.
The shared matrix constitutive owner receives the physical `dt² J M^-1 J^T / 2`
compliance, including cross terms, and owns any native matrix or retained-branch
dispatch. Source coefficients, material identity and state remain with the law.

`RigidTangentialPhase::Dissipative` retains the existing endpoint opposition
and resolved-arrest gates. `StoredEnergy` additionally admits a contact layer
returning stored mechanical energy. The adapter independently checks finite,
nonnegative old/new storage and its difference against the constitutive receipt;
actual body work must still close storage change plus nonnegative exported heat.
Zero storage at both endpoints retains the dry gates even in `StoredEnergy`.
The original directional callback entry preserves its dissipative semantics.

Both routes share one mobility, impulse, momentum/work audit and cancellation
publication owner. Exactly the original pair and accepted history commit
together. Pose/normal evolution and tangent-frame history transport remain
external to this held-geometry phase. `tests/rigid_tangential.rs` covers coupled
full inertia, retained Stribeck history, analytic preloaded-layer energy return,
dry overstep/effort refusals, native-trial cancellation and foreign history.

## Continuous directional friction trajectories on finite rigid bodies

`tangential_friction::directional::trajectory::apply_rigid_directional_trajectory`
composes the original native elliptical relaxation with the actual two-body
mobility. Original masses, full tensors, frozen poses and contact arms determine
the response. Integrated opposing impulses update exactly the selected pair;
the native curve supplies net slip, unsigned travel, physical heat and any
finite arrest time.

Per-body and total linear/angular momentum, actual endpoint velocity and
kinetic work are reconstructed independently. With constant symmetric
mobility, integrated impulse paired with mean endpoint velocity equals kinetic
change even when force turns; this identity is independent of the trajectory's
net displacement. Momentum-increment energy evaluation preserves small losses
under large common translation. Normal recoil is exposed to its original owner.

Both bodies and the original native history publish together after every
physical gate and the final cancellation poll. Other bodies and all poses are
unchanged. The phase holds normal load, geometry and law coefficients fixed,
with no external tangential drive, normal constraint, material query or thermal
store. The native quadrature's estimated-error and unresolved-time limits remain
visible. Three new G1/G3/G4 groups cover full-tensor arrest against independent
impulse-response columns, body-frame covariance, split continuation and atomic
budget/cancellation/geometry refusals.

## Joint Hertz sliding and spin on two finite bodies

`tangential_friction::hertz::apply_rigid_hertz` consumes the original
`HertzWrenchContact`, `RigidContactBinding`, tangent frame, native body momenta,
mass properties and contact arms. It derives all nine entries of the
radius-scaled response `W = J M^-1 J^T`, including translation/spin cross terms;
no effective mass, inertia, prescribed support or supplied slip replaces either
finite body. The original shared Hertz owner receives released displacement
`dt * relative_rate` and held-effort midpoint compliance `dt² W / 2`.

Equal/opposite native point impulses and normal-axis free couples are staged
on only the two bound bodies in fixed-size storage. Actual endpoint point and
angular velocities determine midpoint shear/spin work independently. A fresh
evaluation of the original integrated Hertz law at those actual midpoint rates
checks each force and moment residual and the unpartitioned local heat. Pair
kinetic change, actual work, constitutive work/heat, and total linear plus
orbital/spin angular momentum must close before one atomic publication. The
pair's changes equal whole-system changes; unaffected bodies are not copied or
included in energy subtraction. Cancellation and late geometry/momentum/work
refusal leave the complete original body slice unchanged.

The original law gives feasible global static reactions priority. Zero midpoint
slip does not imply endpoint arrest: a held shear force or spin couple that
would assist either actual endpoint rate refuses as an unresolved slip event.
The trajectory owner must refine that interval; this consumer does not override
static selection, snap endpoint velocity, manufacture an inelastic event, or
turn a residual into heat. Exact kinetic stopping is retained when resolved by
the admitted phase. Normal load/radius, geometry and pose stay held and caller
owned; before/after normal relative speeds remain visible. There is no normal
impulse, external drive, elastic-history replacement, rolling loss, heat store,
contact discovery, time-event search or complete pose integrator in this rung.

Six native G1/G3 groups in `tests/hertz_friction.rs` pass: independent finite-body
translation and spin arrest with support recoil; a manufactured joint Hertz
wrench with full rotated anisotropic cross response and independent native
endpoint work; global rigid covariance, common boost and body-index permutation;
coarse reversal refusal, late angular-momentum failure, cancellation and exact
retry; and static/frictionless limits. The three original rigid tangential and
seven two-axis rolling groups also pass against the same original-source build.
This focused result is separate from a full Cargo/DSR workspace gate.

### Repeated source-coefficient trials on one original rigid pair

`PreparedRigidHertzPair` snapshots only the bound native A/B states, mass
properties, actual contact points/frame, interval, and full radius-scaled W/C.
Repeated coefficient trials reuse those original rates and energy/momentum
totals. Each returns `RigidHertzCandidate` with both native endpoint states and
the original physical report. Normal load and radius remain bit-identical;
original-source static/kinetic coefficients may vary. Source identity and
applicability remain with their original owner. Trials never accumulate another
proposal or change accepted bodies. Actual endpoint work, native Hertz heat and
effort, momentum, reversal and cancellation gates run for every candidate.
`apply_rigid_hertz` delegates to preparation and candidate before its sole final
publication. Three new groups cover independent anisotropic snapshot response,
legacy parity, temperature-coefficient-dependent recoil/heat, held-patch and
effort/reversal/cancellation refusal with exact retry. All three and the seven
prior Hertz groups, including nonprincipal body-coordinate covariance, pass in
the original-source native target.


## Rigid assembly, release and inelastic latching

`assembly::RigidAssembly` combines arbitrary oriented native full COM inertias
and component poses using the parallel-axis law. It computes the composite COM
relative to the supplied reference and caches one admitted full tensor.
Construction and connectivity events are linear in the component and impulse
counts; no dense 6N matrix is assembled. All returned components retain their
original body coordinates and native `RigidBodyState` owners.

Release inherits the parent's rigid velocity field `v_com + omega cross r`
and each component's correct spin before applying prescribed internal impulses
at common physical points. Equal/opposite free couples are separate from point
force moments. Net impulses are accumulated per component before evaluating
work, so multiple impulses sharing a body retain their quadratic cross terms.
Every component survives in the result; spin and orbital momentum are both
included about the same event COM.

Latching is the instantaneous perfectly inelastic rigid-connection law at
unchanged component poses. Total linear and angular momentum determine the
locked velocity. The removed relative kinetic energy is evaluated as a
nonnegative mass-metric squared velocity residual. Accounting uses a common
translating frame, so small locking losses are not hidden by a large common
translational kinetic energy. This loss has no automatic thermal destination.
The caller owns connection geometry, admissibility and the physical latching
decision; no contact discovery, restitution, finite attachment forces,
deformation, fracture, or stored elastic energy model is implied.

Actual returned states supply momentum and energy/work residuals. Declared
part/impulse caps and absolute physical gates apply before publication, and
bounded cancellation leaves borrowed inputs unchanged. Fixed input order is
deterministic; reordered floating sums have only roundoff-level equivalence.
Constant-storage double-double reductions retain linear momentum, rotated spin,
orbital cross-product residuals, and translating-frame kinetic energy until
the signed endpoint balance is formed. Public absolute energies remain rounded
binary64 values and can coincide despite nonzero separation work. A rounded
aggregate must not conceal an actual native endpoint defect: if the single
latched momentum or released component momenta cannot meet the caller's gate,
the event is refused. This does not increase native-state precision or imply
arbitrary-precision conservation.
The rocket staging adapter keeps its dry/propellant partition and inertia
assembly admission and uses this same release kernel with stack workspaces.
It reconstructs receipts after conversion to its stored velocities; the
two-piece rocket path remains allocation-free.

Eleven native assembly cases check an independent full parallel-axis assembly,
zero-impulse rigid motion inheritance, shared-body impulse work against native
sequential impulses, closed-form orbital-to-spin latching, frame covariance,
small loss under large translation, small couple work under large common spin,
native linear/spin defects hidden by rounded totals, cancelling orbital-moment
products, and atomic limits/cancellation/retry.
The existing nine rocket attitude/staging cases remain applicable. These are
floating numerical balances, not exact arithmetic or a complete hybrid-event,
material-calibration, interval-enclosure, or whole-workspace qualification.


## Automatic sustained spherical articulated support

`articulated_contacts::trajectory::sphere_planes::support::advance_supported_carrier_planes`
advances the original balanced spherical carrier through automatically discovered
inelastic landing and sustained support. It reuses the original carrier model
admission, radial material witnesses, simultaneous circular Coulomb response,
and native force/contact/Lie-drift step. Internal actuator effort retains base
recoil; uniform gravity and held link-body torques retain their momentum and
work accounting.

All nearby planes enter one simultaneous solve. Every declared translating
plane, including inactive planes, receives the complete real-arithmetic affine
numerical-model COM-path check plus the actual rounded endpoint check.
Co-located link origins and COMs make the native zero-load momentum-coordinate
update affine independently of joint and attitude evolution. The check uses
the native arithmetic `X0 + R(M*v_body)*(h/M)`. Fractional-time floating
evaluations have no interval enclosure; observed negative gaps still refuse
without projection. A crossing trial halves privately; an accepted approach
within the explicit activation band receives regenerated contact witnesses on
the next step. Plane positions are supplied at the window start and translate
once. No configuration, contact arm, or velocity is projected.

The window exposes accumulated original momentum, external impulses,
actuator/gravity/body-torque impulse work, actuator/gravity displacement work,
prescribed-support work, contact loss, and native drift defects. Physical gates
apply both to individual steps and to the complete window. Late numerical,
geometry, budget, and cancellation failures publish no state. At most 21 planes,
4096 links, 4096 attempted native steps, and one million plane observations are
admitted.

This is first-order inelastic support with finite spatial activation resolution.
Every plane must declare zero restitution; bouncing remains with the original
event owner. Floating affine-model observations are not interval enclosures or
exact event-time certificates. General articulated collision discovery,
arbitrary external link forces, shape-enclosure certification, configuration
projection, and a complete Moreau--Jean method remain outside this owner.
Native contact loss includes supported-load splitting loss and does not by
itself define frictional heat.

Seven focused native cases establish reaction-wheel-driven rolling against
independent constrained mechanics, sliding-to-rolling impulse balance,
automatic landing, discovery of an inactive wall during existing support,
moving-conveyor/elevator work, simultaneous corner support, world-coordinate
and plane-order consistency, and atomic failure/retry.

### Joint native Hertz shear/spin and sphere bulk rolling loss

`tangential_friction::hertz_rolling` resolves five work-conjugate coordinates
`[du_x,du_y,a*dtheta_spin,R*dtheta_roll_x,R*dtheta_roll_y]` with the full native
pair response from actual masses, full inertia tensors, poses and arms.
`PreparedRigidHertzRollingPair` retains only the original pair, source-bound
rolling history, geometry and fixed interval. Its two-dimensional viscous block
is eliminated analytically before the original three-dimensional Hertz solve.
No independent force/spin cap, diagonal effective inertia, fitted damping,
velocity cutoff or numerical compliance is introduced.

The original material coefficient gives `M_roll=-A*N*R*omega_roll`. With
`d=A*N/(R*dt)` and `B=I+d*Crr`, the original Hertz problem receives
`f_h-Chr*d*B^-1*f_r` and `Chh-Chr*d*B^-1*Crh`. Rolling is then recovered and
one combined point/couple impulse updates the two native candidates. The
viscous coupling stays active at zero initial rolling because shear traction
can generate rolling during the interval. Static support, insufficient-traction
sliding, relative spin and both rolling axes share this one physical response.

Actual endpoint momenta independently determine work and heat. Original Hertz
surface loss and original sphere bulk loss remain separate; only bulk loss
enters the rolling checkpoint's cumulative heat, through one original native
advance per accepted phase. No elimination residual becomes heat.
`advance_at_temperature` uses the unchanged source card at actual trial
temperature and preserves its original before/midpoint/endpoint applicability,
frequency, patch and work-coordinate lineage. Repeated admitted coefficient
trials never accumulate another candidate's history or body state.

Positive-load phases refuse initial or induced normal motion because this owner
supplies no normal reaction. Endpoint assisting shear/spin/rolling work requires
interval refinement. The exact released-load and zero-viscosity/friction limits
remain explicit. Source, effort, work, compatibility, momentum, energy and final
cancellation gates precede simultaneous body/history publication. Poses, normal
evolution, partial-slip memory and heat storage remain separate owners; the
existing conservative partial-slip work-overlap refusal is unchanged.

This couples the unchanged symmetric Hertz traction law and original
normal-pressure deformation-loss law under an explicit normal/shear decoupling
approximation. It omits traction/normal-field feedback, including perturbations
of friction integrals by the rolling pressure asymmetry. It is not a complete
first-order mixed-contact continuum solution. Passivity and conservation claims
refer to these two named models inside their retained validity limits.

Seven new native groups cover independent planar static/sliding and induced-roll
references; original Hertz and unconstrained-viscous limits; mixed motion with
independent spatial pressure quadrature; world/gauge/index covariance; endpoint
source/geometry/event admission; and loaded checkpoint, late audit, budget and
cancellation retry. All seven, the ten existing Hertz groups and the seven
existing pure-rolling groups pass against the current original-source
full-inertia implementation. These analytical synthetic fixtures are separate
from measured validation and a full Cargo/DSR workspace gate.

## Partial separation into retained rigid subassemblies

`assembly::partition::AssemblyPartition` prepares a complete partition of an
existing rigid assembly. Every original component belongs to exactly one
nonempty group; group IDs are contiguous, and original component order is
retained within each child. Child mass, COM, complete rotated inertia and
parallel-axis contributions are admitted once by the original assembly owner.
Preparation/storage are O(parts + groups), with no dense multibody matrix.

At release, separation impulse endpoints name original components. Each pair
must cross a released interface and acts on the entire retained group at the
original parent-COM-relative physical point. The original release kernel owns
inherited rigid velocity fields, equal/opposite force impulses and free couples,
full-inertia recoil, shared-group work cross terms and actual output momentum /
energy gates. Events cost O(groups + impulses). There is no intermediate release
of every leaf followed by an artificial inelastic relatch.

Returned group states use the parent's reference axes. Each cached child can
subsequently move independently, be partitioned again for further staging, or
be expanded into original material-body states. Expansion retains child-local
poses, world spin and `v_com + omega cross r`, checks each actual child balance,
and restores original component indexing using O(parts) work and two reusable
largest-group workspaces. Failed admission, balance or cancellation publishes
no partial state and never mutates caller states or the prepared partition.

The six focused partition cases check interleaved full-tensor groups against
original material motion; analytic off-center two-stage recoil and work;
independent native impulse application with shared-group work cross terms;
independent free flight followed by a second separation; single-group and
singleton limits; and malformed membership, caps and late cancellation/retry.
The original seven assembly and nine native rocket rigid-mechanics cases remain
passing in the focused native build. This is not a whole-workspace, target-ISA,
real-time or physical-calibration qualification.

The caller owns the physical connectivity decision, collision/geometric
admissibility, and separation impulses. Retained groups remain perfectly rigid.
This API does not discover fracture, model finite attachment forces or stored
elastic energy, partition heat, or certify the exact continuum trajectory.



## Simultaneous native Hertz contact islands — v1

`tangential_friction::hertz::multiple::PreparedRigidHertzIsland` freezes the
original contact geometry, normal loads, radii and held interval and captures
each participating native body once. Stable unique patch IDs determine
assembly and reduction order. Sparse two-body factors construct the full
`3N × 3N` shear/spin mobility from actual masses, full inertia tensors,
orientations and contact arms, retaining every shared-body cross response.
Untouched spectators neither enter the energy account nor consume the
participating-body budget. Keyed coefficient trials reuse the same original
momenta and prepared response and require exactly unchanged N and radius.

The original constitutive owner now admits dependent positive-semidefinite
contact coordinates without invented compliance. Several distinct patches
between the same two finite bodies can therefore share their translational
and rotational motion. Static load redistribution, unresolved-feasibility
refusals and deterministic nonunique root selection belong to that original
owner; the native adapter does not choose separate per-patch physical states.

Every trial sums all point impulses and free spin couples and forms each body
endpoint once. Every patch evaluates original Hertz traction and heat at its
actual common endpoint-derived midpoint motion. Native force/moment, work/heat,
endpoint compatibility and loaded normal-motion gates apply per patch.
The whole island separately checks actual kinetic change against impulse work
and summed physical dissipation, and checks linear and angular momentum.
Only a complete admitted candidate reaches the publication boundary; failures
and cancellation leave all original body states and preparation unchanged.

This is a held-geometry, held-effort phase. Normal impulse, pose evolution,
retained elastic contact memory and rolling deformation remain with their
original owners. Endpoint-assisted shear/spin and loaded normal motion refuse
unresolved events requiring trajectory refinement. It supplies no thermal
store or experimental validation. Focused G1/G3/G4 groups in
`tests/hertz_multiple.rs` check independent three-body force/spin recoil and
heat, full-tensor response against a separate dense inverse, world/gauge/boost
covariance, exact input-order canonicalization, energetic spectators,
coefficient reuse, single-patch reduction and atomic late-failure retry.
Additional same-pair tests use distinct coaxial and coplanar patches with
actual dependent motion, unequal loads/radii/coefficients, independent summed
impulse and heat, released patches and unchanged cancellation/effort refusal.

## Spatial point loading and finite articulated connections

`articulated::motion::spatial_loads::step_free_floating_spatial` advances the
original free-floating tree under held inertial point forces/free couples and
finite central spring/damper connections between material points on distinct
links. Every point retains its original link frame. Applied world vectors stay
world-fixed while their material application points move; connection length,
axis and strain rate are recomputed at both native Lie-midpoint stages. General
revolute/prismatic/helical trees, noncoincident origins and COMs, and full native
inertia tensors retain the existing ABA and locked-inertia momentum map.

Connections obey `T = k (length-rest_length) + c length_rate`, acting as equal
and opposite central forces with each endpoint's actual moment. They enter the
original generalized acceleration solve but never become an external spatial
impulse. Relative strain rate is assembled from shape velocities and base
rotation directly, so a common translational velocity is not subtracted from
two large endpoint velocities. Both external moments and connection geometry
are assembled before world translation. No dense generalized mass matrix or
new mechanical solver is introduced: work is O(links+loads+connections),
with O(links) additional wrench storage. The original held-wrench API continues
to borrow its wrench slice. Point-only loading skips the endpoint storage pass.

`FreeFloatingMotionStep::internal_body_work` separates connection work from
external and actuator work; it is zero for the original held-wrench owner.
`SpatialMotionStep` additionally reports actual initial/final spring storage,
nonnegative midpoint damping-loss quadrature, connection work/storage defect,
and the complete kinetic-plus-spring work residual. The loss quadrature is not
an exact continuous integral; its thermal destination is unspecified. This
explicit second-order smooth lane requires steps resolving its attachment
modes and does not establish finite-step nonlinear passivity, energy exactness,
rigid constraints, collision response, or a stiff-limit solver.

Admission retains native joint/inertia/state limits, at most 4096 links and
16384 total point loads/connections. Rest lengths, stiffness and damping are
finite and nonnegative. Invalid endpoints, same-link connections, and
coincident/nonfinite attachment geometry at initial, midpoint or returned
states refuse. New passes poll cancellation every 64 entries and the final
poll precedes publication. Caller state and definitions remain unchanged on
late refusal or cancellation.

Seven new native cases independently check a full-inertia spinning body's
world-fixed point-force moment and work; exact damped two-mass modal response
and recoil under timestep refinement; changing-direction connections across a
mixed articulated tree; frame covariance and distant admitted origins; loss
under a 1e9 m/s common axial boost; bit-identical empty spatial loads; and
invalid/late/cancelled refusal with exact retry. All eight existing continuous
articulated trajectory cases also pass. Focused native checks are not a full
Cargo/DSR or cross-ISA qualification.


## Explicit Cayley work chart for native rigid coupling

`RigidBodyIntegrator::step_cayley` uses the same original full-inertia Euler-top
midpoint solve and translation as `step`, while explicitly selecting a Cayley
attitude update for spherical and single-axis rotation too. The default `step`
keeps its original exponential special cases. With the arithmetic rotation
average Rbar, the actual spatial spin impulse is dt Rbar tau_body, and kinetic
work is dt tau_body dot omega_mid. Geometric work owners must pull a spatial
moment back with Rbar inverse and must impose their own finite-rotation chart
limits. No momentum/energy projection or modified inertia is introduced. The
Cayley chart is a second-order angular-phase approximation, not exact spin.
`tests/cayley_work.rs` checks original default behavior, nonprincipal and
spherical inertia, finite kinematic/work identities, spatial torque impulse,
and long free-rotation energy/momentum retention.

## Full-tensor propellant depletion and angular exhaust flux

`rocket::rigid::tensor::{TensorRocketInertia, SpecificInertia, TensorRocketState,
TensorRocketForcing}` extends the original native rocket trajectory to independently
oriented dry, retained-propellant and angular-outflow inertia tensors. The body
frame remains a fixed material frame. `I(m) = I_dry + m K_propellant` can change
principal axes during depletion; no body orientation, momentum, thrust arm or
swirl vector is silently moved into the changing eigensystem. The original
common-axis API retains its exact scalar outflow map and existing default types.

Both models use one generic implementation of the original burnout splitting,
analytic variable-mass translation, torque/swirl kicks and native full-inertia
free rotation. A tensor outflow half-map solves the simultaneous 3-by-3 system
`(I_mid + dm K_exit/2) omega_mid = H_before`, then exports
`dm K_exit omega_mid`. This is implicit midpoint in body angular momentum with
the actual endpoint inertias. It reuses the native fixed-size rotational linear
solve, requires no new dependency or heap workspace, and does not replace a
noncommuting tensor product with three componentwise exponentials. Prepared
specific inertias admit physical semidefinite principal distributions; native
mass properties admit the positive-definite dry and assembled full tensors.
Dry coasting and zero-specific-inertia mass distributions reuse prepared axes.

Returned world angular-impulse/export balances belong to the actual discrete
maps. Smooth forced motion is second order, with timestep refinement required;
implicit midpoint is not an exact noncommuting mass integral or an L-stable
exhaust filter. No complete open-system energy/passivity claim is made. These
are prescribed equivalent thrust, angular co-rotation flux and swirl inputs,
not a derivation of an arbitrary nozzle's internal relative-flow momentum,
Coriolis force, shifting COM, slosh, combustion work or plume dynamics. The
physical distinction between inertia variation and exported angular momentum
follows Eke, *Dynamics of Variable Mass Systems*, NASA/CR-1998-208246,
chapters 2–3: https://ntrs.nasa.gov/citations/19980210404.

Six new focused original-source native test groups and all nine original rigid
rocket groups pass. They cover oblique analytic co-rotation/burnout, independent
noncommuting tensor submap identities, zero-export world-spin retention through
moving principal axes, body-coordinate/world-frame/boost equivalence, and late
budget/cancellation rollback with deterministic retry. A separately integrated
raw-quaternion continuous ODE with a cofactor inertia inverse gives maximum
state errors `2.517403282442654e-5`, `6.292077779579408e-6`, and
`1.5729300554268022e-6` at 8/16/32 steps; its 2048/4096-step references agree
within `2e-13`. This proves the named focused behavior, not a full Cargo/DSR,
experimental validation, cross-platform determinism or performance gate.

## Simultaneous native Hertz and material rolling islands — v0

`tangential_friction::hertz_rolling::multiple::PreparedRigidHertzRollingIsland`
captures unique participating bodies and each original rolling checkpoint,
orders original patch identities canonically, and assembles the complete
`5N × 5N` shared-body response from sparse two-body full-inertia factors.
The original constitutive owner eliminates all material rolling moments
together before the simultaneous Hertz shear/spin solve. Every force and
couple is summed before creating one endpoint per participating body.
Spectators do not enter the island energy or participating-body budget.

Repeated source trials retain exact original load, radius and patch identity
and start from the same captured bodies and histories. The original material
card is admitted at initial and actual midpoint/endpoint rates and its declared
temperature/frequency. Each original rolling history advances once on the
common native endpoint. Original Hertz surface loss and sphere bulk rolling
loss remain distinct. Per-patch force/moment, compatibility, endpoint, work,
heat and cumulative-history gates precede island kinetic-work, total-energy
and linear/angular-momentum closure. Failure or cancellation publishes no body
or checkpoint; the prepared owner remains reusable without hidden warm starts.

One deforming sphere must retain one material parameter/source identity and
held bulk temperature across its contacts. Ordered A/B surface labels may
differ while referring to that same physical sphere. Loaded spherical contact
footprints must be disjoint. Wrapper patch IDs must equal original rolling
patch identities, preventing relabeling of an existing loss history.
Controls bound participating bodies, at most 128 contacts, at most 409600
full-response entries, constitutive sweeps/evaluations and physical efforts.
Dependent full responses use the original PSD Hertz/static-feasibility owner.

Geometry, normal loads and frequency are held. Loaded normal motion and
endpoint-assisted shear/spin/rolling effort refuse unresolved events requiring
trajectory refinement. The parent's first-order normal-pressure rolling loss
and symmetric-Hertz-traction decoupling remains, with independent separated
patch fields. This does not resolve interacting bulk deformation fields,
normal/pose evolution, retained partial-slip storage or experimental validity.

Seven focused G1/G3/G4 groups in `tests/hertz_rolling_multiple.rs` check
independent coupled rolling torque/heat equations, dependent opposite pads on
one pair, friction-induced rolling from rest, original single-patch agreement,
zero/released ports, full-tensor 5N mobility against an independent adjugate
inverse, world/tangent-frame covariance, exact permutation/reuse/spectator
invariance and physical ownership, cancellation, budget and late-failure
refusals. The original seven joint-pair and nine Hertz-island groups also pass.

## Implicit finite attachments on free-floating articulated trees

`articulated::motion::spatial_loads::implicit::step_free_floating_spatial_implicit`
advances the same original `SpatialMotionLoads` with a simultaneous nonlinear
Lie-midpoint solve. It admits the original full link inertias, offset origins
and COMs, mixed fixed/revolute/prismatic/helical trees, held joint efforts,
rotating link-frame wrenches, inertial material-point loads, uniform gravity and
all declared finite Kelvin--Voigt connections. It does not replace the existing
explicit selector or change that selector's arithmetic.

The nonlinear unknowns are only the midpoint joint rates and three base angular
rates. At a solution, `q_h=q_0+h v_h/2`, `v_h=v_0+h a_h/2` and
`R_h=R_0 exp(h omega_h/2)`. The original ABA supplies `a_h`; the original locked
composite inertia reconstructs `omega_h` and base translation from the actual
stage's total spatial momenta. External force and COM-torque resultants depend
only on stage geometry, so their impulse equations are eliminated analytically.
The COM midpoint is `X_h=X_0+h P_h/(2 M)`. Final COM, linear momentum and central
angular momentum use the same midpoint resultants, retaining their world-origin
angular impulse as well as base/shape recoil. Internal connections add no
external momentum. Every connection uses the original material points and
actual stage strain rate; no connection-by-connection velocity update is used.

`fs-solver::nonlinear::NewtonKrylovState` owns iteration and globalization;
`fs-solver::finite_difference::FallibleFiniteDifference` owns the fallible
directional Jacobian actions. The direct `fs-solver` dependency adds no new
package. Native residuals cost O(links+point loads+connections), and bounded
Krylov storage is O(restart*(joints+3)); no dense generalized mass matrix or
Jacobian is formed. The common gravity-only external case skips the extra
external-resultant probe. Numerical Jacobian error is not an adjoint, and its
inner linear solve is deliberately inexact: a 1e-3 forcing floor avoids trying
to solve below directional-difference noise. This does not relax the independently
repeated outer residual or the caller's strict absolute acceptance tolerance.

Controls declare separate physical angular and linear velocity scales; helical
coordinates use rad/s, while prismatic coordinates use m/s. The nonlinear norm
is dimensionless after those scales. A result retains actual unscaled joint and
base-angular equation defects, the repeated scaled residual, native Newton
history and the number of real residual evaluations. Hard limits retain 4096
links, 16384 combined spatial loads, 64 Newton attempts, 64 Krylov vectors, eight
restart cycles per attempt and 16384 residual evaluations, including stencils
and final verification. Cancellation polls all native evaluation and generic
solver boundaries. Malformed controls, failed source/stencil evaluations,
unresolved equations, exhausted work, invalid final geometry and late
cancellation return no advanced state. A forbidden nonlinear trial refuses;
there is no clipping, hidden smoothing, or changed constitutive law.

Storage, nonnegative midpoint damping loss, connection work defect and complete
kinetic-plus-storage work defect use the original shared spatial-load account.
They remain separate from external/actuator work and from any thermal
destination. In the fixed-axis quadratic oscillator limit this is the ordinary
implicit midpoint rule: A-stable, but **not L-stable**. Its scalar stability
function is `(1+z/2)/(1-z/2)` and tends to -1 as a negative real `z` becomes
large. The original MOOSE implementation documents the same linear distinction:
<https://mooseframework.inl.gov/docs/doxygen/moose/classImplicitMidpoint.html>.
Coarse stable steps can retain unresolved high-frequency oscillation and phase
error. General rotating nonlinear trees are not claimed energy-exact, passive,
variational, symplectic or rigidly constrained. This does not implement
generalized-alpha, contact, loop closure, uniqueness certification or an
adjoint. Only initial, trial-midpoint and endpoint connection geometry is
evaluated; no continuous-path admissibility or collision guarantee follows.

Seven new focused native groups pass alongside all seven existing spatial-load
and eight existing articulated-motion groups. Independent tests include stiff
two-mass Cayley response with `h omega=10` over 60 steps, actual base recoil,
undamped storage and damped loss; both three-mass spring-chain modes with a
shared link; exact continuous damped-oscillator refinement; a noncoincident
rotating full-tensor tree under inertial and rotating body loads; frame changes,
rate scaling and a 1e9 m/s common axial boost; a zero-joint force-driven body;
and early/late cancellation, nonlinear and final-evaluation budgets, malformed
controls, final coincident geometry and deterministic retry. The stiff scalar
cases use at most 12 real residual evaluations per step. Continuous oscillator
state errors are `1.07617e-2, 2.68014e-3, 6.69381e-4`; the independent fine explicit
trajectory comparison on the changing-geometry tree gives `5.97471e-3,
1.49312e-3, 3.73255e-4`, with accumulated work defects `1.21004e-2,
3.04096e-3, 7.61237e-4 J`. These are second-order refinement evidence on the named
fixtures and focused original-source native checks, not a full Cargo/DSR gate,
cross-platform qualification, runtime benchmark or experimental calibration.
The broader `frankensim-ext-mbd-core-vqqt` remains open.

### Finite spatial bushings on original material frames

`articulated::motion::spatial_loads::bushing` adds anisotropic translational and
finite-rotation Kelvin--Voigt elements between distinct original link frames.
`evaluate_bushing_forces` exposes the original body-frame wrenches, storage,
internal power and nonnegative damping loss. Translation strain is the B-origin
offset expressed in A minus its declared rest offset. Rotational strain is the
principal SO(3) logarithm of A-from-B, with the desired rest orientation encoded
by the material frames. Nonnegative diagonal coefficients permit inactive axes
without inventing a small stiffness. This does not infer cross-axis constitutive
coupling, rigid six-axis joints, contact, plasticity, or moving rest commands.

The native inverse SO(3) left Jacobian maps relative angular velocity to the
actual rotation-vector rate; its transpose maps conjugate strain force to the
physical couple. The translation law supplies the additional A-side `d cross F`
couple required by the rotating anisotropic frame. Combined forces and moments
are internal and cancel in world coordinates. Relative rates cancel common base
translation and rotation algebraically before evaluating constitutive forces.
No internal bushing rounding is reclassified as an external momentum source.

`step_free_floating_bushings` reevaluates the same law at original explicit Lie
midpoint stages. `step_free_floating_bushings_implicit` uses the existing native
Newton--Krylov midpoint solve, original ABA and locked-inertia reconstruction;
it does not introduce a second dynamics engine or a generalized mass matrix.
Both compose with held joint efforts, world point loads, gravity and finite
axial connections. Storage and loss fields include all active finite elements;
actual kinetic/storage/work defects remain visible. Empty bushing lists retain
the original integrator results exactly. Evaluation is O(links+elements), with
the existing bounds of 16,384 combined loads and 4,096 links, and cancellation polls.

Active rotational elements refuse the principal-log chart boundary at pi minus
one microradian. Initial, midpoint and endpoint checks do not enclose the entire
continuous path or authorize a traversal through that cut. Pure translational
elements do not impose a spurious rotational chart restriction. Implicit linear
stability does not establish phase accuracy, L-stability, nonlinear passivity or
a uniquely selected nonlinear trajectory. Independent finite-difference energy
gradients, spatial moment/virtual-work identities, large common rigid motion,
analytic stiff translational/torsional recoil, combined-tree refinement and
late-cancellation/retry fixtures exercise the actual public consumers.

`bushing::coupled` adds explicitly supplied full 6x6 passive stiffness and
damping in native `[translation, principal rotation vector]` coordinates.
`PassiveBushingMatrix` retains every original SI coefficient and admits exact
symmetry and semidefinite structure after congruence by caller-declared metre
and radian coordinate scales. Negative/unresolved pivots, coupled zero diagonals
and unrepresentable arithmetic refuse; no symmetrization, eigenvalue clipping
or positive loss floor repairs the data. The factor reconstruction residual is
reported in equilibrated coordinates. Exact zero/free directions are admitted;
near-dependent passive forms may refuse without an exact-rank certificate.

`PassiveBushingMatrix::scaled` multiplies an already admitted form by a finite
nonnegative physical coefficient. It retains the original Gram directions,
scales them by the coefficient's square root and rechecks native coefficients,
coordinate congruence and equilibrated reconstruction without refactorization.
Exact free directions stay free; a positive coefficient or factor cannot vanish
silently through underflow. Unit scaling preserves the original representation
exactly and explicit zero scaling supplies the zero form. This permits calibrated
temperature-dependent Kelvin viscosity to reuse the original passive shape with
fixed-size, allocation-free work. It advances no physical state and supplies no
parameter work for a change in elastic stiffness. Seven focused matrix tests,
including three scaling regressions, cover dependent translation/rotation modes,
independent cross work, very weak positive directions and representability refusal.

`CoupledSpatialBushing` fixes those prepared K/C forms to the original material
frames. Full cross-axis and translation/rotation forces use the SAME strain,
log differential, reciprocal wrench lift and finite-state admission as diagonal
elements. Non-diagonal Gram actions supply nonnegative elastic storage and
viscous power; diagonal actions retain the original products exactly. Both
explicit and implicit stepping entry points reuse the original native integrator,
including simultaneous coupling, base recoil, external spatial loads, physical
loss and separate numerical work defects. The diagonal public API is unchanged.
This native law does not infer missing material claims, thermal feedback,
active/circulatory forces, moving rest commands or a different rotational chart.

Four focused preparation tests and three new native integration groups pass
with the original numerical owners. They cover full cross work and physical
rescaling, dependent/free modes, near-active and asymmetric refusal, independent
finite-rotation energy/rate checks, reciprocal force/moment balance, large common
rigid motion, actual coupled trajectories and loss, exact diagonal parity for
both native integrators, and final cancellation. All six original bushing groups
also pass. This focused execution does not constitute a full-workspace gate.

### State-dependent internal effort in implicit motion

`step_free_floating_spatial_implicit_with_effort` and
`bushing::step_free_floating_bushings_implicit_with_effort` evaluate internal
joint effort at the actual reconstructed midpoint in every native Newton
residual and directional-Jacobian probe. Each callback starts with the held
effort vector and may add to or replace it. The accepted midpoint effort also
supplies actuator work. Original ABA and locked-inertia reconstruction retain
carrier recoil; internal effort contributes no external spatial impulse.

Smooth joint springs, dampers and feedback therefore participate in the same
solve as all spatial connections. One O(joints) effort buffer is reused across
evaluations; held-effort calls allocate no additional buffer. Existing work
budgets, physical limits, cancellation, and final residual acceptance remain.
The linear midpoint stability limit still requires timestep refinement for
response accuracy, and arbitrary active effort is not assumed passive.
Constitutive storage and loss interpretation belong to the supplied law.

Callbacks must return the same effort for the same inputs and may not advance
accepted controller/material history during probes. Their side effects cannot
be rolled back. They run at elapsed time h/2 in solver-dependent order, with
cancellation checks before and after each call.

Five focused tests cover stiff reaction-wheel damping and carrier recoil,
inertial-frame speed feedback with analytic second-order refinement,
simultaneous joint/bushing spring-damper work, positive active work, and
callback/domain/limit/cancellation/budget refusals. The 27 existing focused
implicit-connection, bushing, distance-loop and explicit-effort tests also pass.
These selected-source native runs use the original numerical owners; full
Cargo/DSR qualification and the broader multibody bead remain outstanding.

### Simultaneous implicit external forcing and internal actuation

`step_free_floating_spatial_implicit_with_forcing` resolves smooth environmental
forces and joint efforts together with native free-base and joint motion.
Its callback receives the actual trial pose, body twist, joint coordinates and
rates at elapsed time h/2. Reusable link-body wrench and effort slices start
from the held loads on every residual and directional-difference evaluation.
Wrenches remain about each original link origin in its body frame; internal
efforts retain equal/opposite carrier recoil and contribute no external impulse.
Original world point loads and finite axial connections share this trial state.

Three COM world-velocity unknowns augment the original joint/angular solve;
trial geometry and translational drag therefore participate in the same native
matrix-free Newton--Krylov equations. The accepted load sample drives ABA,
momentum and work accounts without another callback. The remaining COM-velocity
residual is returned explicitly. Its change-of-origin contribution to central
angular momentum preserves that sample's actual world moment even at finite
solver tolerance; no endpoint projection or invented force repairs the result.
Held, effort-only and existing constraint paths retain their original n+3
unknowns and evaluation arithmetic.

Linear midpoint stability permits stiff drag but is not L-stability or an
exact relaxation law; nonlinear response still requires timestep refinement.
Callbacks must be smooth and repeatable, and must stage their own history until
acceptance. Native physical limits, finite arithmetic, bounded residual work
and cancellation apply before publication. Constitutive storage, thermal
allocation, discontinuous contact and callback adjoints remain with their own
owners; finite-step work defects remain visible.

Seven focused tests cover independently solved stiff COM drag, positional
force/damping, simultaneous external braking and internal motor recoil,
second-order nonlinear drag trajectories, original held-connection parity,
finite-residual world torque and work, and atomic callback/limit/budget/
cancellation refusals. All seven and five original implicit-effort regressions
pass with the original native dynamics and numerical solver owners.

## Environmental forcing in closed free-floating mechanisms

`articulated::motion::spatial_loads::implicit::distance` now resolves smooth
state-dependent environmental wrenches inside the original distance/point-loop
solve. The diagonal and full-matrix bushing forcing entry points also compose
with constant affine joint rows through `distance::joint`. The original
nonlinear system gains three midpoint COM-velocity coordinates; all material
and joint reaction multipliers remain in the same simultaneous solve.

The forcing callback receives the actual complete trial state at elapsed time
h/2, including carrier pose, velocity and the current initial constraint impulse.
Native loaded ABA, inertia responses and material geometry remain authoritative.
The final verification reapplies the actual summed native constraint impulse
and samples the force law at that reconstructed midpoint. Its retained sample
owns external impulse and work; finalization does not call the law again.
Callbacks must be smooth and repeatable over numerical derivative probes,
including probes outside the constraint manifold, and must stage their own
history until the entire step succeeds.

`ForcedBilateralLoopStep` retains the original loop receipt, accepted complete
midpoint and physical COM-velocity residual. The canonical continuation state
is the original `state_after`, after the final native tangent impulse. Both
constraint events preserve actual finite carrier recoil and signed quadratic
work. External work, genuine actuator work, passive storage/damping and
constraint work remain distinct; a constraint-work defect is not heat.
The original position, tangency, mobility, physical-limit, budget and cancellation
gates remain in force. Held-load and effort-only paths retain their n+3 equations.

This adds smooth environmental feedback to existing rigid closures, not contact
activation, nonsmooth force histories, a nonlinear-branch enclosure or an
energy-exact/symplectic constrained integrator.

Four new physical/refusal cases pass with 28 existing native loop regressions.
Independent reduced equations check coupled carrier translation/spin, full
bushings, environmental drag, actual material/joint reactions and signed work.
A nonlinear closed triangle converges at second order against an independently
assembled full-inertia reduction. Diagonal/full paths, final-sample ownership
and atomic late refusal are also covered. Both changed original test targets
typecheck against the full actual production dependency graph.

## Rigid bilateral distance-loop motion with native reactions

`articulated::motion::spatial_loads::implicit::distance` advances freely moving
articulated trees with declared positive rigid distances between original link
material points. `step_free_floating_distance_loops` admits the same external
loads and finite axial elements as the implicit lane;
`step_free_floating_distance_loops_with_bushings` composes rigid loops and the
finite spatial bushings in one simultaneous solve. Its empty-bushing path is
exactly the plain operation. No mass is clamped, rigid distance is not represented
by a large spring, and a finite bushing's storage or damping is not erased by a
simultaneous rigid loop.

The step uses a SHAKE/RATTLE-style position/hidden-constraint construction:
initial central constraint impulses, the original implicit Lie midpoint, and
final central tangent impulses. Initial impulse magnitudes and midpoint rates
are solved together with the existing Newton--Krylov/numerical-Jacobian owners;
there is no nested nonlinear trajectory integration. Cached native inverse-mass
columns give the initial joint-rate response. Total momenta still use the native
locked-inertia map; internal central constraint impulses have no physical net
external force or couple. The converged initial magnitudes are then applied by
`apply_free_floating_impulses`, and the nonlinear equations are independently
rechecked using that actual summed response and actual spatial momenta. The
endpoint velocity correction uses the original free-floating mobility, including
all link inertias, COM offsets, levers, joints and base recoil. It also passes
through the original impulse owner. Positions are never snapped or projected.

Input positions **and** velocities must already satisfy the requested physical
distance/tangency tolerances. Both are checked again on the final returned
state. The endpoint tangent solve is a bounded dense Cholesky of the constraint
mobility, not a dense generalized mass matrix. The mobility must be reciprocal,
positive and independently reacting: diagonal response is compared with the
same free endpoint bodies, and the unit-diagonal matrix must have resolved
Cholesky pivots. The chosen dimensionless floor must be at least 4096 machine
epsilons. Redundant, already imposed, insufficiently regular or singular loops
refuse with the stage and constraint index; the implementation does not select
arbitrary redundant reactions. These floating admissions are not certified
rank/conditioning or a manifold-distance enclosure.

At most 32 loops and 4096 links are admitted. Initial response preparation costs
O(loops*links), with O(loops*joints) cached storage. Each nonlinear residual adds
one original FK endpoint pass to the native midpoint work. Only the at-most
32-by-32 constraint mobility is factored; endpoint response columns are not
retained. Original bounded Newton work, scaled physical rates, stencil refusals
and cancellation remain authoritative, including the final actual-impulse
verification. Late nonlinear budget exhaustion, invalid final geometry, failed
endpoint regularity or cancellation returns no advanced state.

`DistanceLoopStep.state_after` is the final constrained state. Its `midpoint`
field explicitly contains the intermediate smooth step before the final tangent
impulse. Both full native impulse responses, signed per-loop impulses, actual
position and velocity residuals, nonlinear history, consumed residual work and
mobility pivots are retained. The initial and final impulse work use actual
average pre/post event velocities, including the quadratic impulse term. Their
sum remains visible as `constraint_impulse_work`; it is a finite-step constraint
work defect, not friction heat or an omitted energy source. The complete
`mechanical_work_residual` is actual kinetic-plus-finite-attachment storage change
plus modeled damping loss minus external and actuator work. It closes against
the smooth residual and both native impulse work/residuals. Actual world linear
and orbital-plus-spin angular momentum defects are reported independently.

The construction follows the distinction between primary and hidden constraints
in McLachlan, Modin, Verdier and Wilkins, *Geometric Generalisations of SHAKE and
RATTLE*, <https://arxiv.org/abs/1207.3367>. The underlying native Lie midpoint is
not asserted symplectic, so this composition does not inherit that theorem or an
energy-exactness claim. Intermediate drift states need not lie on the constraint
manifold. There is no continuous-path or collision guarantee, automatic assembly,
nonholonomic rolling, general gear/joint loop, redundancy resolution, nonlinear
branch certificate, interval enclosure, physical calibration or adjoint.

Eight new focused groups pass: an independently solved loaded rigid slider with
both exact reaction impulses; a closed rotating triangle checked against an
independently assembled full-inertia one-joint reduction; nonlinear circular
loop material geometry and velocity reconstructed in world coordinates; two
simultaneous shared-link loops and permutation; scene covariance; inconsistent,
redundant and already imposed constraints; late cancellation, final-evaluation
budget and deterministic retry; and rigid circular support combined with
anisotropic finite bushing stiffness/damping. Triangle state errors refine as
`3.09730e-6, 7.74030e-7, 1.93190e-7`, with mechanical work defects
`2.72326e-8, 7.00372e-9, 1.73603e-9 J`. Forty nonlinear circular steps have maximum
observed world gap `2.22e-16 m` and range rate `2.82e-16 m/s`; their nonzero net
constraint work `6.40818e-6 J` is retained. The combined bushing case retains
`0.0994186 J` damping loss separately from `1.14822e-4 J` mechanical defect.

All 36 combined native cases pass: eight new loop groups, six new spatial-bushing
groups, seven existing implicit-connection groups, seven existing spatial-load
groups and eight original articulated-motion groups. These are focused
original-source checks, not a full Cargo/DSR gate, cross-platform proof,
performance qualification or general multibody lane completion. The broader
`frankensim-ext-mbd-core-vqqt` remains open.

### Implicit actuator effort with rigid bilateral loops

`step_free_floating_bilateral_loops_with_effort` composes admitted distance
and point loops, finite bushings and repeatable state-dependent internal joint
effort in the same native solve. Every residual/Jacobian probe and the
independent actual summed-impulse verification evaluate the law on the
reconstructed midpoint, including the trial initial reaction. The accepted
effort supplies actuator work. Original free-base recoil, both impulse-work
accounts, endpoint closure/physical tangency and mobility admission remain
authoritative.

One effort buffer is reset from held efforts and reused; existing held-effort
paths retain their original evaluator. Trial callbacks must not commit
controller/material history, and arbitrary active effort is not assumed
dissipative. This extends the bounded midpoint lane without a new solver,
coupled electrical state or stronger finite-step energy/uniqueness claim.

Three focused native groups verify driven constrained sliders against an
independent reduced mass of 1.5 kg, including both reaction impulses and
translating-carrier recoil; simultaneous joint-spring/motor and follower-bushing
work; and trial/final-verification faults, cancellation, limits, budget refusal
and deterministic retry. All eight existing distance-loop and seven point-loop
groups also pass. Full Cargo qualification remains unavailable for the
previously reported sibling dependency gap.


## Hertz moving phases and physical collinear stop/restart — v0

Original Hertz pair and simultaneous-island consumers, with and without
material rolling, now supply accepted phase ownership to the original law.
Any nonzero accepted shear/spin retains its kinetic capacity; it cannot spend
a larger static capacity to choose a premature zero-midpoint displacement.
Existing held-effort endpoints still require their original physical audits
and may refuse an unresolved stop. A velocity tolerance is not an instruction
to erase an initially nonzero motion.

`hertz_rolling::collinear::PreparedRigidHertzRollingCollinear` resolves actual
longitudinal slip arrest and subsequent static or rolling-driven restart for
one finite body pair. Preparation retains original full-inertia mechanics and
requires an invariant longitudinal-shear/lateral-rolling subspace, including
zero unresolved normal response. With q=[u,R*omega] and original response
W=[[a,b],[b,c]], gamma=A*N/R, sliding obeys
`F=-mu_k*N*sign(u)`, `r'=b*F-c*gamma*r`, `u'=a*F-b*gamma*r`.
The original exponential path integrates shear travel, force and moment
impulses, surface work and bulk squared-rate work independently.

A bracketed true zero of shear velocity ends a moving phase. Only there is
`F_static=(b/a)*gamma*r` tested against the original static capacity. If it
holds, surface heat is zero while rolling decays with
`r'=-gamma*(c-b*b/a)*r` and continues heating the original bulk material.
Otherwise the outgoing acceleration starts the appropriate kinetic phase,
including reversal. Stop searches retain arithmetic-scaled velocity residuals
as well as caller time/velocity gates; even real 1e-20 m/s initial slip is
integrated rather than declared stationary.

The opaque continuation binds a physically located stop to its exact original
body states, masses, full inertias, body indices, arms and contact frame.
Every phase advances the original material-owned checkpoint; the final body
candidate applies one summed native wrench impulse and independently checks
rates, kinetic work, history heat and total momentum. Preparation/trials are
immutable. Source, effort, phase/root budget, numerical audit or cancellation
failure returns no partial body, stop marker or material history. Work is
bounded by at most 4096 phases and 256 root iterations per event.

This owner holds poses, contact points, normal load, radius, material card and
frequency over the window. It supplies no normal impulse or
normal elastic storage, moving geometry, general shear/spin event trajectory
or simultaneous-island stopping path. Thermal trial owners may refresh admitted
coefficients while reusing the same immutable geometry. Seven G1/G3 groups in
`tests/hertz_rolling_collinear.rs` use independent scalar/exponential references
for true arrest, post-stop rolling, reversal, split continuation, tiny motion
and atomic source/domain/inertia refusals. Two additional `hertz_friction`
regressions verify high-static translation/spin and zero-kinetic free motion;
existing pair and island covariance/recoil/work tests remain active.

Restart extremum bracketing uses the original deterministic logarithm with a
compensated positive `log(1+x)` evaluation, retaining weak departures without
platform libm. The reversal regression continues through the second physical
arrest and subsequent sticking, checking both stopping times against an
independent scalar exponential/bisection reference and continued bulk loss.

`PreparedRigidHertzRollingCollinear::new_driven` retains that same constitutive
path under explicit held world forces and body torques on the original pair.
Native inverse mass and full inertia determine the two driving accelerations;
relative acceleration outside the admitted shear/rolling subspace refuses.
Common translation and rotation retain their own body response and work.
With drive `[d_u,d_r]`, sliding uses `u'=a*F-b*gamma*r+d_u` and
`r'=b*F-c*gamma*r+d_r`. During sticking the actual required reaction is
`F_static=(b*gamma*r-d_u)/a`; the original exponential path locates its first
static-capacity exit and starts kinetic motion in the outgoing direction.
A zero-slip departure marker remains tied to the accepted native state and
the same held drive, capacity and material damping. Changing those inputs
re-evaluates static admission. No momentum is projected to manufacture rest.

Driven stop searches split at the single possible velocity extremum, exposing
an interior stop even when the unresolved interval endpoints have the same
sign. They retain arithmetic-scaled local velocity residuals and the existing
total root-iteration budget. The zero-viscosity limit uses its exact linear
stop, and the original undriven constructor keeps its existing event path.
Full-body external work is integrated through the contact/complementary
mass-metric decomposition. Independent native energy and world momentum
checks subtract that work and the commanded external impulses; neither common
acceleration nor actuator work becomes friction heat. Surface friction and
the original material-owned rolling bulk heat stay separate.

Five driven integration groups add independent constant-acceleration force,
torque and heat references; hidden interior arrest; motor-driven static exit
with split continuation; common-acceleration heat invariance; and tiny-slip
and equal-static/kinetic departure boundaries. They also exercise incompatible
drive and cancellation refusal. These remain held-geometry contact-momentum
trajectories, without gyroscopic or changing-contact geometric acceleration.

### Frictional impacts with closed rigid distance linkages

`articulated_contacts::friction::constrained::resolve_constrained_articulated_frictional_contacts`
resolves supplied Coulomb impacts together with independent positive rigid
distances between original material points, within or across free articulated
trees. It reuses the original articulated impulse response, friction law and
actual endpoint/energy/momentum checks. The empty-distance entry delegates to
the original contact owner. Supplied distance geometry and range velocity must
already satisfy the declared physical tolerances; this event does not locate
collisions, assemble mechanisms or repair positional/velocity drift.

The native augmented impulse mobility is `[D B; B^T A]`. The small normalized
SPD bilateral block A is factored once with the original `fs-la` Cholesky owner.
Cached recovery columns produce the constrained contact mobility
`D - B A^-1 B^T`; the original non-associated Signorini--Coulomb solver supplies
contact impulses j, and the distance reactions are `-A^-1 B^T j`. Contact and
bilateral impulses are accumulated together before the original native update
of each involved tree. Joint/base recoil, full inertia tensors, lever arms,
cross-contact response and cross-tree force transfer all remain in that map.
This follows auxiliary-constraint response condensation, as described in
Baraff, *Linear-Time Dynamics Using Lagrange Multipliers* (1996), section 8,
<https://www.cs.cmu.edu/~baraff/papers/sig96.pdf>.

Positive distance impulses separate the named B point from A; their forces are
equal, opposite and collinear at the actual frozen geometry. Both orbital and
spin angular momentum therefore use the same internal wrench pair. Every
actual returned contact velocity and range rate is checked after the summed
native application. The receipt exposes distance errors/rates, reactions, and
actual constraint work `sum(lambda*(rate_before+rate_after)/2)` in joules.
Admitted input tangency defects remain in that work and the original global
energy residual, which may refuse them; they are not hidden as physical heat.
The moving-support work and contact determinacy policies remain unchanged.
A velocity-only contact witness does not imply unique distance reactions.

The constrained entry admits at most 21 contacts, 32 independent distance rows,
128 trees and 4096 total original links. Vanishing mobility relative to isolated
endpoint bodies, nonreciprocal or insufficiently resolved normalized mobility,
redundant constraints, invalid geometry, exhausted bounds and cancellation
refuse atomically. Only the small bilateral block is factored; the native mass
response stays linear in the involved tree size per row. There is no dense
generalized mass matrix, sequential contact kick or later velocity projection.

`tests/articulated_constrained_impacts.rs` checks analytic unequal-body impulse
sharing and distinct tangential mobility; an independent fixed-joint reduction
with full tensors, COM offsets and rotated scenes; simultaneous shared rods and
redundant contact velocities; smooth rigid-loop motion before and after impact;
exact empty-constraint delegation; and geometric/rank/budget/cancellation
refusal with deterministic retry. The existing contact and contact-step suites
exercise preservation of the original friction and prescribed-support paths.
This supplies instantaneous distance-constrained impacts, not continuous
cross-tree integration, spherical/weld impact closure, CCD, finite-patch spin
friction, compliant history or a unique-reaction guarantee.

### Shared bilateral impact geometry and inertia condensation

`friction::constrained::operators` owns the original objective distance,
spherical and weld covectors and the normalized bilateral Schur construction.
Both native constrained entry points consume this owner; the empty and
all-distance delegation paths remain intact. Borrowed native kinematics admit
fluid consumers without cloning trees.

The mass provider supplies one augmented contact/bilateral mobility and an
uncancelled positive reference response per bilateral row; normalization and
`fs-la` Cholesky preserve the existing independent-row admission. Geometry
retains the A-frame point-closure couple and equal/opposite weld couples, so
admitted nonzero closure error does not create an artificial internal spatial
impulse. The helper supplies no mass, contact law, drift correction, endpoint
authority or unique local reactions on its own.

### Spherical point closures mixed with rigid distance loops

`articulated::motion::spatial_loads::implicit::distance` now admits
`PointCoincidenceLoop` and `BilateralLoop::{Distance, PointCoincidence}` through
`step_free_floating_bilateral_loops` and its `_with_bushings` entry point.
A spherical closure keeps two original material points coincident while
permitting relative rotations. Distance-only APIs delegate to the same engine;
their original config, result and error names alias the common bilateral types.
Mixed closures share the original nonlinear midpoint equations, native floating
mobility, actual summed reaction impulses, finite bushing storage/loss and work.

Point rows use exact objective coordinates `g = R_A^T (x_B-x_A)` and
`gdot = R_A^T (v_B-v_A-omega_A cross (x_B-x_A))`. An impulse lambda in A's
material frame applies `F=R_A*lambda` to B, `-F` to A, and the A-side couple
`-(x_B-x_A) cross F`. That couple is the true differential of the moving
material frame. It enters the same mobility, native impulse application and
quadratic work account, preserving internal linear and orbital-plus-spin
angular momentum even for an admitted nonzero gap. It vanishes at coincidence.
It is not a later angular-momentum correction or a penalty torque.

`BilateralRow` maps each scalar reaction back to its input constraint and point
component. `PointClosureResidual` retains material-frame separation/rate and
actual world separation/relative velocity. Initial and final point norms,
objective-rate norms and physical relative-velocity norms must all satisfy
their explicit tolerances. A small gap spinning rapidly cannot pass by having
zero material-frame rate alone. Initial/final geometric couple impulses and
their actual event work are returned; that work is already included in the
complete reaction work and must not be counted twice.

Cached initial native response columns carry their actual rounded force and
COM-moment resultants into the trial momentum coordinates. Acceptance still
reapplies the accumulated native impulse and reevaluates the actual midpoint
equations. Results expose external-only momentum defects, actual reaction
resultants, and balances after subtracting those resultants. Endpoint impulses
use the endpoint mobility and are checked against the actual returned state.
No geometric projection, dense generalized mass matrix, frozen base or second
trajectory solver is introduced.

The existing bounds are 4096 original links and 32 independent scalar rows;
each point closure consumes three rows. Free-endpoint mobility and normalized
SPD pivot admission refuse already-imposed, redundant, dependent or unresolved
rows without dropping components or assigning arbitrary reactions. Nonlinear
work budgets, geometry/source failures and cancellation remain atomic. This
regular independent-constraint slice does not supply automatic assembly,
general hinge/weld/gear closure, nonholonomic rolling or constraint redundancy.
Intermediate drifts need not stay on the manifold; finite reaction work and
integration defects remain explicit, without general exact-energy or
symplecticity claims.

The seven integration groups in `tests/articulated_point_loops.rs` compare an
independently reduced full-inertia ball mechanism, multiple simultaneous point
closures plus a distance, a nonlinear crank/carriage closure with a rotational
bushing, finite-gap geometric-couple work, covariance and physical-rate/rank/
budget/cancellation refusal. The private row-law unit case independently
differences all twelve endpoint pose directions and checks common-motion
invariance and internal force/moment cancellation. All 44 focused cases passed,
including the 36 existing distance/bushing/spatial/motion regressions.

## Driven joint Hertz shear, spin and rolling

`PreparedRigidHertzRollingPair::new_driven` and
`apply_rigid_hertz_rolling_driven` admit held world forces and body torques in
contact A/B order through the original five-coordinate shear/spin/rolling
owner. Native full-inertia acceleration supplies released drive displacement
before the existing rolling Schur elimination and joint Hertz solve. Original
accepted shear/spin rates still select the kinetic/static phase. Each body
receives one summed external/contact resultant.

Full-body midpoint actuator work and commanded world momentum transfer stay
separate from surface friction heat and material rolling bulk heat. The original
zero-drive entry, source applicability and normal-motion checks remain. An
assisting endpoint effort refuses for refinement. This retains the existing
frozen-geometry held-effort approximation; general multidimensional arrest
events and gyroscopic/pose integration are not added. Four new native test
groups cover mixed powered motion against independent pressure integration,
off-center nonprincipal inertia and frame covariance, exact zero-drive/common
acceleration invariance, and atomic reversal/normal-motion/cancellation refusal.

## Native simultaneous Hertz loading and retained Mindlin shear

`tangential_friction::mindlin::advance_rigid_mindlin` resolves normal indentation
and retained scalar Mindlin shear from two actual rigid-body states, their
masses, full inertias and contact arms. A bounded monotone solve uses the
source-owned exact Hertz mean force and analytic endpoint derivative. The same
affine normal path then enters the original source-bound compliant Mindlin
solve. One combined action/reaction impulse supplies both native endpoints;
no caller-specified effective mass or replacement contact spring is accepted.

Actual midpoint displacement, original normal and shear storage, irreversible
microslip heat, native kinetic energy and world linear/angular momentum are
checked together. Preloaded shear may return energy to the bodies. Source rate
and native shear-speed limits apply at both actual endpoints, so cancellation
in mean motion cannot hide an out-of-domain velocity. The immutable candidate
contains both bodies and complete original-source contact history; any refusal
or cancellation preserves the accepted inputs.

The admitted branch has circular Hertz geometry, a common material point,
central normal force lines and one longitudinal shear coordinate. Transverse
slip, torsional slip, noncentral normal/shear coupling and separation-crossing
intervals refuse. The enclosing owner still supplies initial indentation,
material-point geometry, pose integration and physical-time refinement. G1/G3
checks compare native recoil and work with original source responses and an
independent Hertz ODE, including second-order time refinement. G0/G4 checks
cover unsupported motion, source/history identity, endpoint speed limits,
iteration bounds, cancellation and identical retry.

## Driven native Hertz loading and retained Mindlin shear

`tangential_friction::mindlin::advance_rigid_mindlin_driven` admits held world
COM forces and body-frame torques in original A/B order. Actual mass,
orientation, full inertia and contact arms determine the external point
acceleration. This acceleration enters the released normal and scalar-shear
displacements before the original Hertz mean-force and retained Mindlin solves.
Each body receives one summed external/contact resultant; the zero-drive entry
retains the established action/reaction arithmetic.

Full-body midpoint external work and commanded world momentum remain separate
from reversible normal/shear storage and physical microslip heat. This includes
negative actuator work when a preloaded contact releases energy. The original
source endpoint gates, central geometry and immutable candidate publication
remain. The phase holds geometry and wrenches fixed and assumes affine contact
paths; pose/gyroscopic evolution, separation and recontact remain with their
surrounding owners.

Four new physical regression groups cover clamp-maintained compression and
powered partial slip against the classical virgin Mindlin force and work,
preload retention and actuator energy recovery, nonprincipal full inertia and
world-frame/common-acceleration invariance, exact zero-drive reduction, and
unsupported-motion/source-limit/late-cancellation refusal. All nine native
Mindlin integration tests pass.

## Source-owned adhesive sliding and twisting on finite bodies

`tangential_friction::adhesive` applies the original ordered junction-shear
claims over the actual circular JKR, DMT or Maugis solid-contact disk. Both
finite bodies retain their own mass, orientation, inertia and material arms.
Their complete three-coordinate mobility couples both tangent directions and
normal-axis spin. The uniform junction law shares one force/moment capacity
set, including adhesion at zero or tensile net normal load; its pure torque
is `2*tau*pi*a^3/3`. Neither a compressive-load cutoff nor the Hertz pressure
distribution replaces that adhesive law.

`prepare_rigid_adhesive` stages equal/opposite point impulses and free spin
couples. `apply_rigid_adhesive` publishes both momenta only after independent
native velocity, work, dissipation, energy and total-momentum checks plus the
final cancellation poll. The held wrench must not assist endpoint slip or
spin; an unresolved stopping/reversal event requires a shorter physical
interval. Detached zero-area responses produce exactly zero impulse and heat,
even inside a DMT attractive gap. No cohesive annulus is counted as a solid
shear junction.

The five groups in `tests/adhesive_spinning.rs` use independent pure-spin and
combined disk-rim solutions, all mobility cross terms, finite support recoil,
tangent-frame covariance, stopping-event refusal, detached attraction and
atomic cancellation/audit retry. This held-geometry impulse owner leaves
normal evolution, poses, elastic microslip, peeling and heat partition to
their separate physical owners; the returned junction loss is unpartitioned.

`tangential_friction::adhesive::multiple` extends the same source-owned law to
simultaneous contacts on shared native bodies. Canonical physical patch IDs
fix assembly/reduction order. Sparse inverse-mass Gram factors retain complete
inertias, actual COM-relative arms and every cross-contact term; dependent
wrench coordinates remain PSD rather than receiving invented compliance.
One summed impulse updates each participating body once. The Hertz and adhesive
owners share these native mechanics and work audits, while their original
constitutive surfaces remain distinct.

`PreparedRigidAdhesiveIsland` retains every original normal/shear source and
real disk, including detached sources. Attached patches require coincident
points and closed normal motion at both native endpoints, including at zero or
tensile normal load. Detached zero-area sources create no constraints, wrench
rows or heat; detached-only bodies are excluded from the active energy sum.
An all-detached group returns an exact empty mechanical no-op plus its original
detached identities. Spectators never dilute microscopic work/momentum gates.

The actual final midpoint motion is independently evaluated through each
original uniform-junction integral. Local force, torque, native recoil, work
and heat must agree before the common original-body energy/momentum audit.
Endpoint assisting work refuses an unresolved stop/reversal interval. Prepared
state and source receipts remain immutable through trials; application publishes
every body only after the final cancellation gate. Sweeps are numerical work,
not physical time or heat. Bounded feasibility/convergence limits of the shared
constitutive solver remain explicit; trajectory refinement belongs to the caller.

Four G1/G3/G4 groups in `tests/adhesive_multiple.rs` use two distinct symmetric
disks with independent exact rim-slip integrals, combined point/free-couple
torque, finite support recoil and complete relative-energy loss. They also
exercise source/order preservation, dependent mobility, energetic spectators,
detached/all-detached groups, held-normal admission and atomic event, momentum,
cancellation and budget refusals. `tests/hertz_multiple.rs` remains the original
native mechanics regression target for the shared extraction.
The symmetric microscopic-pad fixture uses a resolved short physical interval
at unchanged force/work gates. Its much larger interval explicitly exercises
support-budget refusal; strong translation/spin conditioning is not hidden by
relaxed tolerances. A separate axial disk isolates the physical stopping-event
refusal without that numerical conditioning limit.

`PreparedRigidAdhesiveIsland::advance_thermal` reuses the captured Gram factors,
full compliance, body endpoints and source-owned disks with the junction claims
refreshed by `SourceThermalAdhesivePatch`. Exactly one source per attached patch
is required; canonical IDs, original normal/shear source identities and exact
disk geometry are checked before a solve. Detached patches cannot acquire rows
or heat through this API. The original fixed-source `advance` and thermal path
share the same numerical solve and native audits, with actual endpoint traction
reevaluated through the selected temperature-dependent uniform-junction law.
Only shear strengths change here: normal thermodynamics, poses and heat transport
remain with their own owners. The additional `adhesive_multiple` regression
checks identical constant-claim mechanics at different temperatures, source
ordering, foreign-card/changed-area refusal and retry from unchanged endpoints.

`tangential_friction::adhesive::trajectory` accepts the original source-owned
fixed-ray relaxation through a finite stop and any remaining zero-drive rest.
It prepares the actual complete pair mobility and applies only the analytically
integrated point/free-couple impulses. Body velocities are never overwritten or
clipped. Actual native endpoints must agree with the located event; per-body and
global momentum, original uniform-disk force/moment, moving-phase heat, normal
closure and full kinetic work are independently audited. Stationary dwell does
not multiply the moving traction power by the full requested duration.

Both `prepare_rigid_adhesive_trajectory` and its atomic application retain the
original source geometry and strengths. Detached zero-area sources remain exact
no-ops even for separated moving points. Poses, normal evolution, external drive,
thermal changes and curved velocity paths remain outside this held-geometry
branch. Four focused groups in `tests/adhesive_trajectory.rs` cover finite support
recoil through pure-spin arrest at zero/tensile load, a full-mobility coupled rim
ray, tangent-gauge covariance, gap detachment, curved-path refusal, cancellation,
late native angular-momentum refusal and deterministic retry.

## Powered shared-body Hertz shear, spin and rolling assemblies

`PreparedRigidHertzRollingIsland::new_driven` admits held world COM forces and
body-frame torques in complete original body order; nonparticipating entries
must be zero. Each body's acceleration is captured once and supplies every
attached patch's released five-coordinate displacement before the original
full cross-contact mobility, rolling Schur elimination and surface solve.
All patch impulses and that body's external resultant are summed before one
native state update.

Independent midpoint external work and commanded world momentum transfer count
each shared body once and remain separate from surface and material rolling
heat. Original source identities, patch separation, material-temperature
agreement, native normal-motion and endpoint-assisting-work gates, immutable
trials and exact zero-drive behavior remain. Geometry and normal loads are
held; general multidimensional arrest, pose evolution and gyroscopic motion
retain their surrounding owners.

Four new native groups check independently solved coupled rolling equations
including dependent two-pad support, nonprincipal full-inertia mixed motion and
common acceleration, original body indexing/patch permutation/singleton limits,
and atomic invalid-drive/normal/slip-event/cancellation refusal. All eleven
native multi-contact Hertz/rolling tests pass.

## Native adhesion after prescribed normal thermal evolution

`tangential_friction::adhesive::thermal_normal` consumes the original JKR
held-indentation thermal refresh, accepted branch token and actual source-owned
disk. It refreshes the same selected junction-strength claims and reconstructs
the complete native response at the new solid-contact radius. The new radius
therefore changes both torque capacity and the scaled mobility; an old prepared
matrix is never reused across changing normal geometry. Zero/tensile resultant
loads retain their actual attached area.

The native solver independently checks body impulses, momenta, kinetic work and
physical shear/spin loss. The combined storage check includes the actual change
in normal elastic-plus-surface potential and the signed parameter-work supply:
delta(K + U_normal) + friction loss - parameter work = 0 within declared bounds.
Parameter work is never exported as friction heat. The atomic application moves
both native momenta, new normal/shear source owners and branch token together
after final cancellation admission, without cloning the new normal owner.

Normal approach, body poses and mass properties are held during this prescribed
thermal then impulse split. Normal-motion integration, expansion, caloric/entropy
closure, changing attachment and held-load actuator motion are not supplied.
Three G1/G3/G4 groups in `tests/adhesive_thermal_normal.rs` check independent
zero-approach JKR radius and torque, separate bulk/interface temperatures, full
scaled mobility, complete storage work, zero-slip thermal work, accepted source
continuation and atomic history/control/cancellation refusal.

### Native curved adhesive trajectories

`tangential_friction::adhesive::transient` applies the original source-owned
curved disk trajectory to the actual finite body pair. It retains every
translation/spin mobility entry and applies only integrated original point and
free-couple impulses. Actual endpoint velocity/spin, each body's requested
impulse, total momentum and native kinetic work are independently checked.
Original instantaneous power quadrature and the physical finite-extinction
tail supply heat; an average-wrench endpoint-power test is not substituted for
the actual changing traction history. The returned stopping interval and
uncomputed travel bounds remain available to the surrounding time owner.

A conservative original normal-response-row bound covers the entire moving
path and analytic tail. Unsupported normal/shear coupling refuses without
inventing a normal correction. Detached zero-area sources remain exact no-ops.
Preparation is immutable and application publishes both bodies together only
after all original dimensional gates and cancellation checks. Source properties,
normal geometry and poses remain frozen with no external drive. Three native
groups in `tests/adhesive_transient.rs` exercise full-mobility curved recoil,
actual per-body impulses and work, tangent covariance, eventual rest with finite
common motion, source gap detachment and atomic normal-path/cancellation/budget
refusal and retry.

## Shared-body assemblies with retained Mindlin partial-slip histories

`tangential_friction::mindlin::multiple::PreparedRigidMindlinIsland` resolves
multiple original source-bound scalar Mindlin histories on shared native bodies.
Preparation canonicalizes unique contact IDs while retaining original A/B
material order and whole-array body indices. Actual masses, full inertia tensors,
orientations and contact arms determine every cross-site response. The complete
physical two-tangent Gram mobility remains available for inspection. Its actual
longitudinal subspace enters the existing bounded simultaneous tangential owner;
the second bookkeeping coordinate is exactly inactive, without regularization.
Patches sharing a native body must have attachment separation greater than their
summed original footprint radii and geometry tolerance. This conservative gate
refuses overlapping or duplicate patches without inventing a contact locator.

`new_driven` includes held world COM forces and body-frame torques in released
contact motion. Every body receives one summed contact/external resultant, and
its independent midpoint external work and commanded world momentum count once.
Each original constitutive evaluation starts from the same accepted checkpoint;
nonlinear sweeps do not advance time, wear or heat repeatedly. Contact storage
and microslip/gross-sliding dissipation retain their original physical owners.

Actual native endpoints must satisfy the scalar phase's normal, lateral and
torsional motion gates. Actual midpoint displacement must match each solved
longitudinal increment and leave excluded coordinates within their dimensional
limits. Per-contact work/storage/heat, complete native kinetic energy and world
linear/angular momentum close before all body and history candidates return
together. Foreign source history, unsupported physical recoil, exhausted budgets
or late cancellation leave the accepted assembly unchanged.

The phase holds original normal loads, material conditions and geometry fixed.
Normal prestress/equilibration, independent patch ownership, pose integration
and history transport remain with the caller; no normal impulse/work, rolling
memory or inter-patch bulk elastic field is inferred. G1/G3/G4 tests compare a
three-body/two-contact load line with independent mass/inertia algebra, classical
virgin and reversal Mindlin work, retained partial/gross/reversal histories,
products of inertia, common acceleration and world/local frame covariance,
canonical contact/global body indexing, and physical recoil/transactional refusal.

## Native DMT normal motion with actual junction sliding

`tangential_friction::adhesive::normal_motion::PreparedRigidDmtPair` consumes the
original source-bound DMT normal/shear load line on two finite native bodies.
An explicit fixed approach datum and actual collinear geometry markers
authenticate the original normal state. Paired forces act at the markers'
common midpoint reaction plane, preserving action/reaction and the work-dual
normal/shear Jacobian. That plane is the declared contact discretization, not a
solution for unequal-material surface-deformation partition. Actual masses,
full inertias and reaction arms supply every entry of the three-point-force
mobility; no caller effective mass or omitted normal/shear coupling is accepted.

The original source integrates reversible normal work and real solid area along
the same affine discrete path, including attachment and detachment. Finite-gap
attraction can exchange native kinetic energy and normal storage while adding
exactly zero shear area. Actual midpoint body motion is replayed through the
original source before acceptance. Differential native kinetic work, original
normal storage, physical junction heat, individual impulses and complete pair
linear/angular momentum have independent dimensional gates. Static and kinetic
strengths stay distinct. Unsupported solid-contact torsion and unresolved
sliding capture/reversal refuse without clipping motion or inventing heat;
gap-only spin remains free.

Preparation and retries leave accepted bodies/history unchanged. Application
publishes both momenta and the source-authenticated endpoint occupancy together
after final cancellation admission. The fixed-Jacobian phase returns actual
COM, geometry-marker and reaction-point midpoint travel but does not advance
poses. The surrounding pose/contact owner must consume that travel and refresh
the geometry before continuation; an advanced normal state against stale marker
positions refuses. The affine constitutive path supplies neither a validated
continuous trajectory nor contact detection, finite rotations or an impact law.
Six G1/G3/G4 groups in `tests/adhesive_normal_motion.rs` cover analytic finite-gap
recoil/storage and pose-owned continuation, integrated-area attachment and
detachment with all cross responses, static support, tangent gauge/body indices,
and source/geometry/torsion/sliding/cancellation/budget refusal with retry.
Native execution remains subject to the existing dependency build limitation.

## Native JKR normal attachment and pull-off

`tangential_friction::adhesive::jkr_normal_motion::PreparedRigidJkrPair`
condenses actual masses, full inertia tensors and reaction arms into the
original `fs-contact` JKR normal trajectory. Geometry markers and an immutable
approach datum authenticate accepted source history. Equal/opposite forces act
at the common midpoint reaction plane. The original impulse and its first time
moment restore native momenta and COM, angular and marker travel through smooth
motion, capture, reversal and pull-off. A whole-interval midpoint impulse does
not replace that event-resolved travel.

`advance_with_wrenches` and the matching prepare/apply APIs accept the actual
held COM force and body-frame COM torque for each original body. Their full
inertial response supplies the scalar trajectory's held load, while tangential
and rotational motion remains unconstrained. Contact and external impulses are
summed before each native momentum is rounded once. The external work is force
times actual COM travel plus torque times actual angular travel; it includes
the contact first time moment and the external acceleration's half-step-squared
travel. Endpoint midpoint velocities do not replace that work. The original
zero-load APIs use the same execution path with zero wrenches.

Independent checks close native kinetic work, integrated original contact work,
normal storage plus explicit snap release minus external work, individual
combined impulses and total pair linear/angular momentum against the actual
external wrenches at the held geometry. Snap release supplies no impulse
and is not automatically assigned to heat. Successful preparation is stack-only;
original scalar substep/root/event budgets bound the trajectory. Cancellation or
late refusal publishes no body or history update; application publishes both
body momenta and the source branch together.

This remains frictionless JKR with one fixed normal Jacobian. The contact supplies
no tangential force or free spin couple. The pose owner must consume travel and refresh
geometry before continuation; stale geometry refuses. No finite-rotation pose
integration, contact discovery, hidden support, junction shear or measured
specimen calibration is claimed. The six G1/G3/G4 groups in
`tests/adhesive_jkr_normal_motion.rs` exercise exact detached flight, analytic
initial attraction and snap/rebound energy, full off-axis inertia, first-moment
travel, world-frame/body-index covariance, geometry-owned continuation and
source/budget/cancellation refusal with deterministic retry.
Four loading groups in `tests/adhesive_jkr_normal_motion/loading.rs` add exact
zero-load API equivalence, common gravity, tensile preload and original pull-off,
off-axis full-inertia force/torque response, first-moment external work, body-index
covariance, and atomic load/domain/cancellation refusal followed by valid retry.

`initialize_motion` and `continue_motion[_with_wrenches]` retain the original
sealed source/mass/policy/time/step token across native phases. New phase loads
may change; actual refreshed geometry and normal rate must agree with accepted
history, and changed effective mass or numerical policy refuses. No clock reset
or event replay makes a mismatched continuation pass. Two continuation groups
compare a split native attachment/rebound exactly with the original scalar
trajectory and reject unconsumed geometry, altered velocity and changed policy.

## Prepared thermal trials of adhesive transients

`PreparedRigidAdhesiveTransient` captures the original pair, point kinematics,
actual radius-scaled mobility and complete normal-response row once. Each
`advance` starts from the same captured endpoints. `advance_thermal` requires
equality of the complete original source-owned patch and uses only its pinned
temperature-refreshed junction law; normal source, disk geometry and body poses
stay held. The original trajectory, integrated impulses, native momentum/work/
energy and cancellation checks share one execution path with the existing
prepare/apply APIs. The normal-path bound uses the refreshed strength and actual
stopping interval. Detached patches retain exact no-op behavior and cannot
acquire a thermal solid face. Three groups in
`tests/adhesive_transient_thermal.rs` cover independent finite-inertia thermal
arrest/heat, retries and original-API equivalence, full source/geometry refusal,
cancellation/detachment, and strength-dependent normal-path/budget refusal.

## Native full-stick mixed-mode JKR normal/shear motion

`tangential_friction::adhesive::jkr_mixed_mode::PreparedRigidJkrMixedModePair`
embeds the original coupled full-stick JKR potential and the original ordered
static junction strength into actual paired-body motion. The complete three
point-force Gram matrix retains normal/shear cross response from both masses,
full inertia tensors and actual reaction arms. The source owner solves all three
coordinates together using its analytic discrete-path tangent, and bounds the
entire constant-effort quadratic phase against source, radius and static-slip
limits. Shear changes real area and normal reaction; no held old disk or diagonal
effective-mass approximation replaces that coupling. A prepared pair reuses the
original small LU workspace across trials.

Normal approach is authenticated by collinear geometry markers and an immutable
normal datum. The two elastic shear coordinates are retained constitutive
history, driven by actual velocities at a common midpoint reaction plane. They
are not tangential gaps between separated markers. All three point reactions
therefore cancel internal force and angular impulse. The native history also
pins the original tangent frame, normal datum and material marker arms; the
source token retains its source, static capacity, full mobility, policy, time
and step count. A changed frame, marker association or mobility requires an
explicit history-transfer owner, which is not supplied by this phase.

Held actual COM forces and body-frame torques drive the same coupled solve.
Contact and external impulses are accumulated before each momentum is rounded
once. Actual midpoint COM/angular travel closes external work for this one
constant-effort phase. Independent native checks reconcile all three rates and
coordinate increments, individual impulses, total linear/angular momentum,
contact work, and kinetic-plus-source storage minus external work. Elastic
unloading may assist body motion and return stored energy; it produces no
automatic friction heat. Cancellation or any source, domain, nonlinear or
native refusal publishes neither body momentum nor elastic history.

This is a fixed-Jacobian, held-pose, full-stick phase. The caller owns pose
integration and any later contact/frame migration. Relative torsion is screened
at both endpoints, which bound its affine rate within this phase. Torsion,
sliding, radius instability and changed mobility refuse without invented spin
traction, detached landing state, snap energy or heat. The five groups in
`tests/adhesive_jkr_mixed_mode.rs` cover the original zero-shear normal limit,
elastic unloading and area recovery, independent full-tensor three-axis recoil,
external force/torque work, world/body-index covariance, held mixed preload and
chronology, and atomic frame/source/mobility/geometry/twist/load/cancel refusals
with retry.

`PreparedRigidJkrMixedModePair::continue_at_temperature` authenticates the
original held-indentation thermal receipt and its pinned shear-strength
successor while preserving the captured body poses, momenta, full mobility,
reaction arms, tangent frame, normal datum, elastic shifts, actual rates,
controls, time and step count. Native history and torsion checks precede source
transfer. The replacement pair and retained token return together with the
original source's separate normal-elastic, tangential-elastic and surface
energy changes; their complete signed parameter work must also pass the native
energy tolerance. This endpoint parameter change applies no impulse, heat or
elapsed time and claims no intervening thermal path, expansion or reattachment.
Three native temperature groups cover unequal source temperatures with nonzero
actual momentum and subsequent whole-body energy closure, exact no-op transfer
after an accepted step, and stale frame/rates/geometry/mobility/twist plus
cancellation refusal with retry.

`continue_at_geometry` accepts compatible caller-owned pose endpoints without
restarting the source trajectory. The private native token seals ordered body
indices, full mass properties, accepted endpoint momenta, accumulated pending
COM travel and whether each pending phase had angular motion. The same material
markers and normal datum must reproduce the retained approach; actual COM
travel, individual linear and world-spin momenta, total orbital-plus-spin angular
momentum and full kinetic energy must pass the original native tolerances.
An attitude or common spin unsupported by any accepted angular motion refuses.
The existing tangent isometry re-expresses elastic shear, and the source owner
retains its source, static capacity, normal coordinate, policy, time and steps
while admitting the new actual full mobility and rates. Transport storage
roundoff remains in the native energy gate; a change of contact-projected
kinetic energy is diagnostic only, with no work or heat attribution.

This handoff is not an attitude integrator or a general moving-contact path.
Noncollinear original material markers, unowned COM corrections, unsupported
torsion and failed conservation gates refuse without projection or relocation.
The caller still owns finite-rotation evolution and path refinement. Four
geometry groups cover two actual loaded normal phases with retained shear,
changed midpoint arms and full-tensor mobility, tangent/world covariance,
accumulated common travel, zero-motion attitude/common-spin refusal and atomic
geometry, momentum, property and cancellation failures with retry.
Source/native execution remains subject to the documented dependency build
limitation.

## Native Hertz--Mindlin dynamics with two-surface fretting wear

`mindlin::wear::advance_rigid_mindlin_wear` and its driven counterpart compose
the original native changing-normal/shear phase with source-bound Archard
inventories. Actual masses, full inertia, contact arms and held external wrenches
determine both trajectories once. Positive mu uses that phase's physical shear
dissipation divided by mu; explicit zero friction uses the original sourced
Hertz mean normal force times actual accepted shear travel. No endpoint load,
offered motion or reversible normal work substitutes for local sliding.

Original body candidates, retained normal/shear history, normal potential,
shear storage/heat, external work and both cumulative wear volumes return
together. A late wear refusal or cancellation publishes none of them. No new
constitutive solve, source refresh, history copy, mesh or heat partition is
introduced. The existing central, scalar, affine-path and compressed-contact
limits remain; volumes diagnose small removal without changing geometry,
mass, material conditions or the native dissipative work port.

Three G1/G3/G4 groups verify exact unchanged native mechanics through elastic
release and driven loading/reversal, separate ordered k/H inventories, the
nonlinear frictionless normal trajectory, original selected-role identity,
overflow after successful mechanics, zero-drive equality and final cancellation.

## Simultaneous Hertz loading and retained shear across shared bodies

`mindlin::multiple::loading::PreparedRigidHertzMindlinIsland::{new,new_driven}`
advances all original normal and scalar-shear histories on a shared native body
island. Actual mass, full inertia, orientations and arms assemble the complete
three-coordinate contact mobility. The original Newton--Krylov owner solves the
full normal load line using exact source-owned Hertz path means and analytic
endpoint derivatives; Krylov products reuse those derivatives without another
source query. The original simultaneous tangential owner then consumes those
same affine normal paths with every longitudinal cross-contact term retained.

Normal and tangential coordinates keep separate physical owners. The supported
mechanical blocks must be orthogonal up to arithmetic roundoff; genuine
normal/shear cross-coupling refuses instead of being discarded. All patches
remain compressed, and original source rate/indentation limits apply to the path
and actual native endpoints. Independent footprint admission covers the largest
radius reached on each affine path, including growth-induced overlap.

Each body receives one summed normal, shear and external resultant. Actual
midpoint displacements, normal potential work, retained shear storage, physical
microslip heat, actuator work, native kinetic energy and world momentum close
before the complete body/history candidate returns. Cancellation and source or
iteration/evaluation-budget refusal preserve all accepted inputs. Separation,
geometry/pose evolution, rotating history and interacting elastic patch fields
retain their surrounding owners.

Four G1/G3/G4 groups verify independently expanded three-body/two-contact load
lines through loading/unloading/reversal, zero-imposed-shear unloading loss,
second-order 8/16/32-step refinement against an independent coupled Hertz ODE,
full-inertia/frame/common-acceleration invariance, canonical body/contact indices,
growing footprints, genuine cross-block and late unmodeled recoil, source/budget
refusal and cancellation. The finest ODE indentation error is below 1e-10 m.

## Spatial forces during variable-mass rigid motion

`RigidRocketState::advance_spatial` applies additional `PreparedSpatialLoads`
through the original common-axis and full-tensor depletion map. A world-held
force at a rotating material attachment uses `r_body cross R^-1 F_world` at
both actual native torque-kick attitudes. Held body and world free couples
retain their distinct frames. The prepared resultant enters the original
analytic changing-mass translation exactly once; callers must not duplicate
these forces or moments in the separate propulsion/external forcing inputs.

Thrust, rotational exhaust, burnout, cancellation and substep work limits retain
their original owners. Prepared load evaluation uses constant storage/work per
substep independent of attachment count. Changing drag, gravity or other
state-dependent laws still require outer-stage evaluation and refinement.
The original `advance` entry selects zero additional spatial loads. Focused
regressions in `rocket::rigid::tests` compare rotating-arm momentum against an
independent planar RK4 solution and check world-couple impulse during mass
and rotational-exhaust depletion (G0/G1).

### Physical exhaust momentum and energy in free translation

`rocket::exhaust` advances physical particle outflow on the original native
variable-mass trajectory. Its comoving mass, outward momentum and complete
energy moments are held during emission; exact burnout and subsequent coast
remain owned by `RocketState`. The inertial exports integrate the changing
vehicle velocity, including both the momentum cross term and advected kinetic
energy. A physical particle momentum is distinct from effective nozzle thrust;
pressure thrust, gravity, external work and rotation are excluded here.

The original thermodynamic owner must debit the complete comoving energy once.
The kinetic frame correction is retained separately from its possibly large or
signed formation datum. Actual endpoint velocity and rounded propellant loss
enter stable momentum and energy residuals. A common boost that erases recoil,
or amplifies a mass-rounding defect, therefore remains visible to the caller's
physical acceptance budgets. Mean velocity and mean squared speed support
accepted-stream rounding corrections without a second dynamics solve. The
mechanical step allocates no heap storage.

G1/G3/G4 fixtures cover independent time quadrature of exported moments,
Galilean transformations, opposed streams, signed energy references, burnout
splitting, tiny depletion, representable extreme velocity variance and actual
inventory/velocity rounding defects. Rustfmt and source review passed; native
execution remains unverified in the checkout missing constellation siblings.

## Point-parcel capture on the native rigid body

`assembly::capture::capture_point_mass` admits one incoming point parcel with
explicit nonnegative mass, body attachment arm and world inertial velocity.
The parcel sticks to the original rigid body without an external impulse.
The result uses the existing `RigidBodyState` and `MassProperties`: body axes
and orientation are retained, the COM moves by `m r/(M+m)`, and central inertia
becomes `I + mu (|r|² 1 - r r^T)`, with `mu=M m/(M+m)`. The returned COM shift
must be subtracted from original body points before native continuation. No
positive parcel radius or fictitious intrinsic inertia is introduced.

The original inverse-inertia Gram factor and bounded 3-by-3 solve determine
the internal capture impulse from contact slip, including full rotational
recoil. Nonnegative sticking loss is evaluated as positive translational and
rotational impulse quadratics and checked against half the slip/impulse work.
Its heat destination remains explicit. Native momentum residuals include both
spin and COM relocation; the energy gate uses a common translating frame and
polarized spin changes, without subtracting large laboratory kinetic energies.
Absolute tolerances remain caller-owned floating-point gates, not certificates.

Positive capture consumes two part slots and one internal impulse slot under
`AssemblyControls`; zero mass is an exact one-part identity. Unrepresentable
mass addition, inertia, momenta, unresolved mobility/work, balance failures and
late cancellation return no endpoint and leave both original owners unchanged.
The event uses fixed-size impulse algebra and the existing bounded central-inertia
decomposition. Six focused test groups
cover analytic central/off-center capture, full-tensor conservation and frame
changes, large common translation, inherited rotating-point capture after
native flight, identity and atomic refusal. This is instantaneous rigid
sticking, without contact discovery, finite-volume flow, slosh, continuous
mass flux, thermodynamic mixing or a complete fluid/body energy model.

### Momentum retained across native motion and load assembly

Free-flight momentum reconstruction and the spatial-load assembly shared by
free and prescribed motion use compensated linear/angular reductions. Held link wrenches, world
point loads, axial connections, and finite bushings retain their low components
until the complete per-link force/couple reaches native dynamics. Shape momentum
remains compensated while combining with or subtracting from base momentum.
Point moments retain both cross-product terms; attachment transport moments and
free couples remain separate through assembly.

Explicit free-flight and prescribed kinetic-work increments use shared link
COM/central-spin polarization. Explicit free-flight gravity work uses the computed
displacement before adding the world origin. These are floating-point accuracy improvements, not exact conservation
claims for unrepresentable state increments or arbitrary timesteps. The existing
integrator, constitutive laws, admission, and cancellation boundaries remain in
force. `articulated_load_reduction` exercises actual supported/free trajectories,
off-center moments, and held/world/axial/bushing cancellation.

## Point-parcel release with retained rigid recoil

`assembly::release::release_point_mass` removes a caller-identified point parcel
already included in the original body mass and inertia. For original mass T,
parcel mass m, retained mass M=T-m and original COM-relative arm r, the retained
COM shifts by `-m*r/M` and its central inertia becomes
`I - (m*T/M)(|r|² 1-r r^T)`. The original native full-tensor admission rejects
nonpositive or physically inconsistent retained inertia; no eigenvalue clipping,
diagonal replacement or fictitious parcel radius is introduced. Admissible
moments alone cannot establish the caller's geometric inventory.

Both outputs inherit the original rigid velocity field before a declared
internal impulse acts on the parcel and its exact negative on the retained
body at the same physical point. The parcel has no intrinsic spin or free
couple. Rotational recoil uses the point arm measured from the retained COM.
Separation work is a nonnegative sum of translational and full-inertia impulse
quadratics; the pre-event point velocities coincide, so initial-spin cross
terms cancel. Its energy source must be debited by the caller.

The retained body uses the original native momentum state; the free parcel
stores canonical world momentum with a checked velocity query. The returned
COM shift must be subtracted from original body points before continuation.
Actual returned momenta and both positions enter linear and spin-plus-orbital
angular residuals, including any rounded first-moment defect. FMA-centered
momenta and polarized spin changes give the common-frame energy/work audit.
These are caller-budgeted floating balances, not exact-state certificates.

Positive release consumes two part slots, and only a nonzero impulse consumes
an impulse slot. Zero mass with zero impulse is an exact identity; zero mass
cannot receive an impulse. Unrepresented depletion, inadmissible inertia,
nonfinite arithmetic, balance refusal and late cancellation publish no result.
The event uses fixed-size impulse algebra and the existing bounded native
inertia decomposition. Six focused G0/G1/G3/G4 test groups cover analytic
central/spinning recoil, full-tensor frame covariance, capture/release
interoperability, large boosts, lost geometry and atomic refusal. This is
instantaneous material detachment; no fracture, contact discovery, finite
parcel extent, slosh, continuous mass flux or thermodynamic source is inferred.

### Physical exhaust under gravity and applied forces

`RocketState::advance_with_exhaust_forced` composes the same physical particle
moments with held `ExhaustEnvironment` gravity and applied force. The original
native changing-mass trajectory and burnout event remain authoritative. Exported
momentum and energy integrate the velocity during emission; the external force
continues through coast. Gravity impulse integrates retained mass, and gravity
work integrates the correlated mass and velocity. Work is mechanical input to
the retained vehicle, with no invented propellant caloric model or potential
energy datum.

The analytic moments use `d=F*b/M0-w*x`, avoiding `F/q` when the expelled rate
is tiny. Positive small-depletion series preserve mean drift and velocity
variance when force acceleration opposes gravity. No numerical quadrature,
heap allocation or additional time steps enter this mechanics path. Zero loads
select the original free-exhaust calculation; a closed or exhausted vehicle can
still move under the declared loads.

Both external impulse and external work are returned explicitly. Residuals
subtract them in a chart about the initial velocity while using the actual
native endpoint, exposing velocity and inventory rounding. Constant loads and
relative exhaust moments are analytic assumptions, not a claim of exactness for
changing thermochemistry, drag, rotation or shifting centre of mass. Independent
time-quadrature regressions cover noncollinear burn/coast, gravity/mass
correlation, force cancellation, small depletion and Galilean work transforms.

## Source-owned directional friction on supported mechanisms

`articulated_contacts::mixed::resolve_mixed_articulated_elliptic_contacts_with`
extends the original source-history transaction to freely moving, grounded and
externally driven articulated trees. It samples each active original source
history exactly once, in canonical contact-ID order, using actual native
material-point locations, normalized contact axes and incoming relative
velocity. Inactive geometry does not query or advance history. Returned source
candidates retain their original input indices; accepted histories are borrowed
and need not implement `Clone`.

The source callback supplies the held preimpact friction ellipse and a private
history candidate. The existing mixed inverse-mass/contact solve determines all
finite-body recoil, coupled joint response, native support reactions and signed
boundary work. The original total-work range witness remains authoritative for
redundant events; a source callback does not turn representative contact loads
or individual source-work shares into identified values.

All mechanical and source candidates publish together only after the native
endpoint momentum, work, speed-limit and constitutive gates pass. A later source
error, mechanical refusal or cancellation returns neither set of candidates and
leaves all borrowed accepted owners unchanged. Callback-owned external side
effects require the caller's own staging. The free-only entry retains its
original default response type and shares the same transaction implementation.

This is a frozen event with coefficients held at the native incoming query. It
does not add a nonlinear velocity-dependent source iteration, source texture
transport, pressure inference, thermal solve or contact-time discovery.

Four focused groups in `tests/articulated_mixed_sources.rs` cover independently
derived full-tensor hinge recoil and driven work under incoming-speed source
queries; canonical history mapping with inactive and redundant driven contacts;
atomic source, speed-limit and late-cancellation refusal; and exact free-entry
parity. All four groups and the five existing mixed-normal groups pass in the
bounded native numerical target. Unrelated tangential, trajectory, bilateral
and executor paths are excluded; blocked dense GEMM is guarded by an unexercised
fail-fast fixture. This does not substitute for the full Cargo or DSR gate.


## Full passive bushings on prescribed supports

`articulated::motion::spatial_loads::bushing::coupled` exposes supported
explicit stepping with and without stage-dependent joint effort. Fixed and
driven mechanisms now use the original complete 6x6 passive stiffness and
damping forms, including translation/rotation coupling. The existing finite
material-frame geometry, SO(3) pullback and reciprocal wrench lift remain the
only constitutive owners.

Diagonal and full-matrix supported adapters share one midpoint lifecycle. The
same forces enter native ABA joint motion and RNEA support recovery, preserving
linear/angular impulse, signed support work, actuator work, attachment storage
and damping loss. No effective mass or substitute fixed base is introduced.
This explicit lane must resolve the fastest coupled mode; sampled chart checks
do not certify the continuous path.

`tests/articulated_prescribed_coupled.rs` checks an independently solved
translation/rotation mechanism under stationary and translating supports,
feedback evaluated at the actual stages, and exact diagonal-law parity.


## Implicit spatial motion on prescribed supports

`articulated::motion::spatial_loads::prescribed::implicit` advances fixed
and driven open-chain mechanisms with the same native implicit-midpoint
solver family used by free flight. The unknowns are original scalar joint
rates. Native ABA evaluates acceleration at each trial's midpoint joint
configuration and caller-supplied base boundary; world loads, axial connections,
diagonal or full passive bushings, and smooth feedback effort share that stage.

The original Newton--Krylov, FGMRES and fallible directional-difference owners
provide bounded iterations, cancellation and residual-evaluation budgets.
No dense mass/Jacobian matrix or surrogate floating base is introduced.
The final residual is independently evaluated before publication. All-fixed
trees still recover actual reactions and work without fictitious unknowns.
Feedback must be smooth and deterministic: nonlinear trials and derivative
stencils may revisit states, and callback side effects remain caller-owned.

Explicit and implicit lanes share admission, native RNEA support recovery,
endpoint continuation, and physical momentum/work accounts. Support impulse
and signed support work remain external; bushing storage and damping loss
retain their original owners. Numerical momentum and energy defects remain
visible. The caller supplies mutually consistent base pose, twist and spatial
acceleration samples. Midpoint is A-stable on a linear oscillator but is not
L-stable, does not select a unique nonlinear branch, and does not certify
continuous-path geometry or rigid closure.

Focused tests include stiff axial and cross-coupled translation/rotation
oracles, second-order resolved oscillation, jointly implicit feedback and
external loads, translating/accelerating/rotating supports, scalar root
actuation, zero-DOF trees, exact diagonal parity and atomic refusal.


## Rigid payload capture on prescribed supports

`articulated::motion::inventory::supported::capture_payload` connects the
original payload inventory to fixed and driven articulated mechanisms. It
combines the incoming full rotated inertia with the target link, maps the
actual incoming-versus-attached momentum mismatch through the original
prescribed-base impulse owner, and returns the loaded model with the same
native supported motion state used by explicit and implicit continuation.
Time, base boundary, joint coordinates and topology are retained; joint rates
respond simultaneously, including scalar root joints.

The support supplies the external linear/angular impulse and may perform
signed work. Whole-system momentum includes both the incoming payload and
retained bodies. Actual per-body increments and a separate zero-bias RNEA
projection audit world and joint momentum. Capture loss is the nonnegative
sum of relative-velocity kinetic quadratics, not the lab-frame energy decrease.
A common inertial accounting frame retains support work and exposes the
remaining numerical energy defect. Caller-supplied absolute residual and loss
budgets, original endpoint speed limits and cancellation precede publication.

The attachment token retains the exact original link inventory. This event
models instantaneous perfectly inelastic attachment at supplied geometry;
contact discovery, finite grasp dynamics, extra gear/loop constraints and a
thermal destination remain separate owners. Independent tests cover fixed and
moving support catches, mobile roots, branched response, full off-centre inertia,
world-frame and inertial-boost covariance, loaded gravity continuation and
atomic refusal.

## Simultaneous material detachment from one rigid snapshot

`assembly::release::release_point_masses` releases several caller-identified
point parcels already included in the same original rigid mass distribution.
All positions and inherited velocities refer to the original snapshot. The
retained COM and full central inertia account for the complete removed mass
and first/second moments, with one original native inertia admission. The
original point-release descriptors and canonical free-parcel type are reused.
Zero-mass entries require zero impulse and produce no parcel; positive outputs
retain their input order. Empty or all-zero input is a validated identity.

Each parcel impulse has an equal/opposite action on the retained body at its
original physical point. Shared recoil is coupled: the retained translation
uses the SUM of all impulses and its rotation uses the SUM of their moments
about the retained COM. Separation work includes the squared sums and their
cross terms, plus each parcel's own impulse energy. Independent single-event
work cannot be added to represent this simultaneous event. Initial-spin terms
cancel because both sides share each original material-point velocity.

Compensated fixed-order reductions take O(inputs) work, one retained-inertia
admission and output-sized parcel storage, with no dense parcel-by-parcel
system. The part budget includes the retained body and EVERY input slot,
including zero masses, before any scan or allocation; nonzero impulses consume
impulse slots. Canonical parcel momenta use FMA, and the retained momentum is
recovered from the actual original momentum minus all actual parcel momenta.
Actual returned positions/momenta, the first-moment defect under translation,
and common-frame/polarized energy enter the existing balance gate. Cancellation
is polled on each bounded pass and before atomic publication.

Seven G0/G1/G3/G4 regression groups cover analytic opposing/co-directed jets,
asymmetric COM/recoil, inherited rigid motion, singleton interoperability,
full-tensor frame covariance and recapture, parcel subdivision/permutation,
large common boosts, lost geometry, budgets and atomic refusal. These are
finite-precision event checks under caller tolerances, without a time-step,
fracture law, geometry-inventory proof, finite parcel extent, heat destination
or continuous thermodynamic stream model.

## Simultaneous point-parcel capture

`assembly::capture::capture_point_masses` attaches many zero-intrinsic-inertia
parcels from one original body geometry and velocity snapshot. It reuses
`PointMassCapture` and `PointMassCaptureResult`; the returned internal impulse is
the aggregate impulse on the original body, and the returned loss covers the
entire event. Every input arm uses the original COM. Subtract the single
returned COM shift from existing body points before native continuation.

Four bounded linear scans use fixed-order compensated mass, first-moment,
inertia, and momentum reductions. Inertia is assembled from positive
parallel-axis contributions about the final COM, avoiding subtraction of two
large global tensors. The original full native inertia owner is admitted once;
all-central input retains its prepared principal frame, and a singleton uses
the original optimized capture path. Batch reduction has constant auxiliary
storage and no dense parcel mobility matrix. Native tensor admission retains
its own existing bounded allocation; no throughput benchmark is claimed.

The final rigid field shares translation and rotation across every parcel.
Its sticking loss is the positive mass-metric residual of the original body's
translation and central spin plus every parcel's velocity. Independently
adding single-hit losses against the unchanged original body would omit this
shared recoil. The endpoint audit reconstructs actual native linear momentum,
spin plus orbital momentum and common-frame kinetic energy, with the actual
COM displacement, first-moment defect and polarized spin changes.

The part budget counts the original body and every descriptor, including
zero mass, before scanning; each positive parcel consumes an impulse slot.
Empty/all-zero input is a validated identity. Individually sub-ULP masses may
be captured when their compensated aggregate is represented. Invalid inputs,
unrepresented aggregate addition, failed balances and bounded cancellation
leave borrowed state unchanged. Six G0/G1/G3/G4 regression groups cover shared
translation, analytic COM/full-inertia/spin changes, corotation, covariance and
permutation, correctly rebased sequential equivalence, singleton receipt
identity, large boosts, lost geometry, aggregate mass and atomic refusal.
Physical inventory, contact/event discovery, finite parcel extent, continuous
mass flow and the destination of removed mechanical energy remain caller-owned.

## Continuous bidirectional particle momentum on the native translating body

`RocketState::advance_with_particle_exchange` admits simultaneous held incoming
and outgoing `ParticleStreamMoments` on the original native state. The signed
mass rate changes its transferable inventory, while both gross streams retain
their original relative momentum and complete energy. Equal incoming and
outgoing mass rates can therefore accelerate the body, exchange energy and
thermalize incoming relative kinetic energy without changing retained mass.

The method integrates `M*v'=F+Pi-Po+M*g` over the complete declared interval.
Inventory below the dry-mass floor refuses; no intake shutdown, outflow switch
or post-exhaustion policy is inferred. The original native loss integrals and
positive velocity-variance calculation also cover gain by time reversal, with
the gain mean coefficient `A-C`. Zero net mass has its analytic constant-mass
limit. Small-fraction series preserve opposing force/gravity increments and
work, and normalization avoids division by a vanishing signed mass rate.

Each `ParticleStreamIntegral` retains gross mass, relative momentum/energy,
inertial momentum/energy, and the separate inertial kinetic correction. Both
streams use the same velocity mean and second moment. External force and
gravity return their complete impulse and correlated mechanical work. Native
residuals use the actual endpoint velocity and separately stored dry and
transferable inventories in an initial-velocity chart. Their signs are retained
change plus outgoing minus incoming minus external input; gross-versus-net
inventory rounding and unrepresented endpoint recoil are exposed.

`retained_relative_energy_j` is the caloric credit required of the composing
material owner: incoming minus outgoing complete relative energy. Incoming
relative kinetic energy is already included once. No second caloric state,
pressure-flow work, source qualification or species budget is invented here.
Held relative intake is an explicitly supplied source moving with the body;
a stationary ambient distribution has changing relative moments and needs a
different source closure. Rotation, COM shifts and finite sticking events stay
with their existing owners, including `assembly::capture` and `release`.

Six focused groups in `rocket/exhaust/exchange/tests.rs` cover independent
physical-flux quadrature for gain/loss/throughflow, constant-mass force and
thermalization, one-way parity and split continuation, closed/zero/cancelled
intervals and depletion refusal, Galilean/datum transforms, and tiny-fraction
and endpoint-rounding limits. Source files were parsed with Rustfmt; these new
native regressions have not been executed in this dependency-limited session.

### Inertial intake with simultaneous external loads

`RocketState::advance_with_inertial_capture_forced` extends the existing
inertial-intake owner with held force and uniform gravity. It integrates
`M'=q`, `(M*v)'=q*u+F+M*g`, the capture position, and all mechanical exchanges
over one accepted interval. Orbital influx follows the moving capture point;
external angular impulse includes the continuously changing retained mass.
Separate force/gravity work and nonnegative inelastic heat accompany the native
endpoint. Heat uses positive relative-velocity moments, never a repaired energy
residual; incoming internal energy and capture heat remain explicit material
credits. Source inventory and incoming inertial velocity remain source-owned.

Zero loads preserve the force-free calculation, zero rate gives ballistic
forced motion, and cancellation/refusal publishes no state. Small-rate series
retain near-balanced force/gravity response; the dominant-intake chart retains
the original body's finite momentum and heat. G1/G3/G4 tests independently
integrate the conservation laws and cover exact endpoint momentum, work/angular
balances, interval composition, frame changes, limiting cases and rejection.
Spin, finite-size attachment, outflow and time-varying loads remain outside this
point-body operator.


## Passive release from supported payload inventory

`articulated::motion::inventory::supported::release::release_payload`
completes supported pickup, native continuous motion, and drop. It restores the
exact retained link inertia from the original attachment token and reconstructs
the released payload's current pose and rigid velocity field from native FK.
The supported time, base boundary, coordinates and rates are copied exactly;
passive detachment supplies no impulse, support work or mechanical loss.

Free and supported release share the original loaded-inventory identity check.
Changed or repeated tokens refuse, and captures nested on one link release in
reverse order. The original full payload tensor is retained; no inertia
subtraction or artificial recoil is used. Only the changed target contributes
to the independent world and joint momentum decomposition audits. Common-frame
energy defects remain visible and must satisfy caller-supplied absolute
allowances. Original state admission and cancellation precede publication.

This removes an ideal rigid attachment without an ejection impulse. Ongoing
contact, finite separation dynamics and an ejection energy source remain
separate models. Focused regressions specify a capture through actual implicit
gravity-loaded carry into release, driven rotational offset velocity, nested
and repeated inventory, and atomic refusal.


## Full passive bushing feedback in free flight

`articulated::motion::spatial_loads::bushing::coupled` now exposes
`step_free_floating_bushings_with_effort` and
`step_free_floating_bushings_implicit_with_effort`. The original complete
passive 6x6 stiffness/damping law and smooth internal actuator feedback use the
same actual native stage, including finite carrier recoil. Original spatial
loads, axial attachments, bushing storage/loss and actuator work retain their
existing owners.

Diagonal and full-matrix explicit entries share one lifecycle; the implicit
entry uses the original bounded nonlinear solve and final residual check.
No additional dynamics pass, controller approximation or constitutive law is
introduced. Trial feedback retains native effort limits, repeatability,
cancellation and immutable-input semantics. Explicit stability and implicit
midpoint phase/geometric limits remain unchanged.

Three focused regressions independently solve a coaxial translating/spinning
mechanism with reduced inertias 6/5 and 10/7, verify both carrier recoils and
complete momentum/work accounts, solve stiff full K/C plus feedback through
an independent 2x2 midpoint map, and retain exact diagonal parity for both
native integrators.

## Automatic sphere encounters across free and prescribed articulated trees

The mixed trajectory spheres owner, through
advance_mixed_sphere_contacts, locates its own sphere-to-sphere event times
and material witnesses, then continues the original free or prescribed
articulated motion. ContactSphere descriptors follow input-system order;
one SpherePairLaw is required for every unordered stable-ID pair. Complete
canonical pair coverage is bounded by 21 pairs, with at most 4096 native links,
64 events and explicit attempted-step/geometry-query budgets.

Each exposed root-attached sphere encloses its original tree. All link origins
and COMs coincide, and joints are fixed or revolute; full physical inertia
tensors and movable prescribed roots retain their native response. Second
moments must fit the declared radius, but do not certify an unknown shape's
enclosure. Free centres have uniform-gravity acceleration with no other link
force. Prescribed centres have constant world acceleration, while their support
frames have constant world angular velocity. Exact rotation samples use the
original exponential owner; body acceleration remains
R^T*a_world - omega_body cross v_body. Angular joint motion, body torques and
relative actuator effort retain the original ABA and RNEA owners.

Every relative centre path is quadratic in real arithmetic. The existing
fs-contact accelerated-sphere owner finds first entry and every interior
distance extremum, including an encounter hidden between separated endpoints.
All actual native endpoint and regenerated root-link material-point gaps must
also be nonnegative. Up to 64 earlier floating-time neighbours are tried;
positions, contact arms, inertia and velocities are never projected. All pairs
within the explicit positive activation resolution enter one simultaneous
native mixed Coulomb event. Discovery restarts from the actual accepted native
velocities and support clock.

Generated witnesses pass the same original contact-frame/normal-line admission
before an event time is accepted. Its unchanged 256-epsilon rounding coefficient
now retains the actual pose-translation and material-arm magnitudes used in
each point transform, so cancellation at the world origin cannot collapse the
arithmetic allowance. Material points and normals remain untouched; meaningful
off-line separation still refuses. The direct geometry target
tests/articulated_contact_geometry.rs covers this scene-frame covariance.

The result retains native event responses, full support reactions, signed
boundary/actuator work, unallocated impact loss and complete reusable states.
Actual stage/event energy and spatial-momentum changes and independent source
impulses are reduced in DD before window balances; large rounded endpoint totals
remain diagnostics. Smooth integration defects remain separate from physical
work and impact loss. Cancellation or a late path, law, model, limit or budget
refusal publishes no partial window.

This is floating analytic coverage of the declared spheres, not a validated
trajectory tube or interval CCD certificate. General offset links, protruding
geometry, arbitrary support commands, sustained contact, Zeno continuation and
impact-to-heat allocation require their own physical owners. Incompatible
prescribed/prescribed approach and unresolved near-rest inward contact refuse.
The focused native target is tests/articulated_mixed_spheres.rs.

## Loaded sustained contact with mixed articulated boundaries

step::mixed::step_mixed_articulated_frictional_contacts composes the original
held-load impulse, one simultaneous native mixed Coulomb event, and unloaded
native free/prescribed drift. It reuses the mixed trajectory state/system
carriers, a common physical clock and explicit midpoint/endpoint prescribed
support samples. A driven support supplies its real reaction impulse and signed
work during the load kick, contact event and unloaded drift; no artificial
support mass or independent link inertia replaces the original coupled model.

The mandatory geometry callback receives the actual post-contact states,
native endpoints and supplied support paths. It must admit the complete path,
omitted-feature coverage and continued applicability of the frozen contact
witnesses. The midpoint and endpoint samples alone do not establish that claim.
Strict original endpoint gaps, joint limits, cancellation and absolute native
work/momentum gates precede publication. All borrowed states remain unchanged
on failure, and callback-owned history must be staged until success.

Complete native load/contact/drift receipts accompany the accepted state.
DD actual stage increments and source reductions retain body/gravity/actuator
kick work, load-support work, drift-support work and independent impulse
balances. Joint and gravity displacement work and their splitting differences
remain separate; prescribed drift kinetic change minus support work is the
numerical drift defect. Contact loss may include removal of a free load kick
and therefore is not automatically continuum friction heat.

Gravity displacement work retains base translation, original relative-FK link
origins and rotated local COM arms as separate DD differences. A distant world
origin therefore cannot erase represented local motion and mislabel its real
work as an additional splitting defect.

The first-order force/contact/drift composition is bounded by 128 systems,
4096 links and 21 contacts. It does not claim a complete Moreau--Jean method,
automatic collision discovery, continuous rotating-body-wrench work, a
continuous angular gravity integral along an arbitrary driven path, or a
thermal partition. The focused native target is
tests/articulated_mixed_contact_step.rs.

### Directional and source-resolved mixed contact steps

`step::mixed::elliptic` uses that same native step for oriented Coulomb ellipses
and original source histories. The circular path/result types retain their
default coefficient and contact-result types. Geometry checks read only the
original witnesses; no alternate load, inertia, drift or work law is introduced.

`step_mixed_articulated_elliptic_contacts_with` queries each active source once
at the actual post-load-kick point velocity and assembled contact frame. Its
non-Clone history candidates remain private through native contact and drift,
endpoint separation, complete-path admission, physical balances and final
cancellation. A late refusal returns no mechanical endpoint or source candidate.
All plane/base work and force/couple reactions retain their native scope.

Coefficients are held for the original first-order discretization. Source
conditions, normal law, texture transport, support paths and geometric coverage
remain explicit; total contact/splitting loss is not a thermal partition.
Focused tests add circular-limit parity, oriented-traction reaction and
8/16/32-step motion refinement, actual post-kick source selection and driven
work, plus source/path/cancellation rollback and retry. The new tests were not
executed in the constrained editing environment.

## Supported finite-payload separation impulses

`articulated::motion::inventory::supported::release::release_payload_with_impulse`
completes the supported rigid-payload lifecycle with simultaneous equal/opposite
point impulses and free couples. It shares the existing `SeparationImpulse`
convention and physical-point transform with free release. The retained tree
uses native prescribed-base ABA and the released payload uses its original full
spatial inverse inertia; every retained joint responds in the same event.

The prescribed support retains its pose/twist and exposes its independent
reaction impulse and signed work. Actual separation work is the impulse pair
paired with the resulting relative velocity jump, admitted against an explicit
source-work budget. Whole-system linear/angular momentum, native joint momentum,
actual kinetic increments and common-frame energy defects remain independent
accounts. Numerical decomposition error is never assigned to loss or source
work. The composing actuator must debit the returned supplied work from its own
energy source.

Passive release and the zero-impulse limit retain identical state and balance
accounts. Native endpoint speed limits remain active; impulse reactions do not
pass through force/torque effort limits. A nonzero impulse whose outgoing payload
jump is entirely erased refuses atomically. This
bounded event does not infer collision, fracture, finite-duration separation or
additional gear/loop constraints.

Six focused cases independently check grounded and moving launchers (including
negative support work), finite launcher recoil, a full-tensor offset payload on
a branched tree, world-frame/inertial-boost covariance, exact passive parity,
source-work refusal, cancellation and actual speed-limit refusal. The full
focused inventory run passes 24 cases: these six plus four supported passive
release, eight supported capture, three free capture and three free release.

## Constant gear and lock constraints in supported implicit motion

`articulated::motion::spatial_loads::prescribed::implicit::constraints`
connects the original affine `JointConstraint` rows to supported implicit
midpoint motion. World and axial loads, diagonal/full passive 6x6 bushings,
and smooth actuator feedback share one actual stage. Each row is constant
`A q = b`, with tangent condition `A qdot = 0`; the native prescribed-base ABA
acceleration receives the ideal reaction from
`(A M^-1 A^T) lambda = -A a_free`.

The original prepared constraint owner now shares its normalized row-space
factorization between free and prescribed boundaries. Prescribed trial columns
use the raw native inverse-mass response, so numerical unit probes are not
mistaken for physical velocity-limit violations. The bounded row mobility is
the only dense matrix; no generalized mass matrix or substitute dynamics is
introduced. Caller feedback and spatial laws are evaluated once per residual,
before the native reaction, with genuine effort limits still enforced.

Original RNEA recovers support force and torque from the completed acceleration.
The common supported motion receipt now carries separate signed
`constraint_work`, zero for unconstrained paths. It is subtracted in the physical
work balance independently of actuator work, support work, passive storage and
damping loss. An explicit absolute reaction-work budget admits any finite
nonzero tangency defect; such work is never hidden in heat or caller actuation.
Initial, accepted midpoint and endpoint affine/tangent residuals must all meet
the original row tolerances. Numerical derivative probes may lie off the
manifold; no clipping or state projection is applied. Empty rows retain exact
original implicit behavior.

Eight new cases check independently reduced geared-rotor inertia and signed
rotating-support work, a coupled translating/rotating rack with full K/C and
feedback, one feedback call per residual, exact diagonal/full parity, a scalar
root lock whose reaction exceeds the motor rating, large virtual response
columns with small real speed limits, row-scaling equivalence, separate finite
reaction work, rank/configuration/tangency refusal, physical endpoint limits
and cancellation. All eight and 39 affected original regressions pass. This
owner covers constant affine gears/locks; activation impulses, backlash,
varying ratios, nonlinear point loops and general nonholonomic constraints
remain separate work. Original midpoint phase and geometric limits remain.

## Free implicit motion with existing gears and locks

`articulated::motion::spatial_loads::implicit::constraints` connects the
original constant affine joint rows to free implicit midpoint motion, including
world and axial loads, diagonal/full passive 6x6 bushings, and smooth actuator
feedback. The native free inverse-mass responses retain the finite carrier
inertia. A bounded row-space solve supplies ideal reactions; the original
momentum-coordinate equations reconstruct the actual carrier recoil and never
add internal reactions to external spatial impulse.

Each nonlinear trial evaluates the caller force laws and loaded ABA once.
Native inverse-mass probes form the row mobility, using the original normalized
prepared factorization. Numerical derivative probes may leave the affine
manifold, while native configuration limits and independently checked initial,
accepted midpoint and endpoint row residuals remain authoritative. No
constraint projection or clipped residual is substituted for the equations.
Empty rows reproduce the original implicit result exactly.

The shared free-motion receipt now separates signed `constraint_work` from
genuine `actuator_work`, connection work/storage, damping loss and external
work. The implicit finalizer uses the existing compensated product reduction
and double-double work arithmetic for caller efforts, reactions, gravity and
body work before rounding the public accounts. Existing distance/ball-loop
entries retain their original callback signatures and solver ownership.

Eight independent cases check geared-rotor reduced inertia 110/9, carrier
recoil and external world force/torque; a coaxial translation/spin gear with
reduced inertia 13/7, full K/C coupling, feedback and both carrier recoils;
diagonal/full equivalence; real motor versus passive-reaction admission;
row-scaling equivalence; finite signed reaction work; degeneracy and
configuration/rate refusal; and cancellation. The native distance, spherical
point-loop and bilateral-feedback regressions also exercise the shared
implicit owner. These are existing constant affine gear/lock relations;
activation, backlash, changing ratios, and joint composition with nonlinear
point-loop constraints remain separate models.

Focused native validation passes 72 cases across 12 integration targets,
including ten new physical/refusal cases and 62 original regressions.

`constraints::step_free_floating_spatial_with_forcing` and the diagonal/full-
bushing counterparts jointly resolve state-dependent environmental wrenches,
caller actuator effort and existing constant joint rows. The original n+6
midpoint owner evaluates the complete trial pose and velocity, calls the forcing
law once, and performs loaded ABA once. The original native inverse-mass row
solve then adds internal reactions; accepted samples own external impulse/work,
caller actuator work and separately signed constraint work. Carrier translation
and rotation recoil and the finite COM residual retain their original momentum
and change-of-origin treatment. Initial, midpoint and endpoint row tolerances
and independent reaction-work admission remain authoritative; no projection
repairs a residual.

Held/effort-only n+3 paths retain their original equations; empty rows exactly
reproduce the unconstrained forcing result, and diagonal/full representations
retain exact parity. Six additional physical/refusal cases pass with 20 original
regressions on the native numerical sources: independently reduced geared-rotor
and translation/spin models with material-point drag, nonlinear locked-rotor
second-order refinement, finite signed reaction work, frame/row rescaling and
atomic refusal. The extended original integration target also typechecks against the
full native dependency graph. These are smooth constant affine constraints;
activation, backlash, varying ratios and new contact/support laws remain outside
this extension.

## Supported payload capture with existing gears and locks

`articulated::motion::inventory::supported::constraints::capture_payload`
connects finite rigid-payload capture to the original constant affine joint
rows. It uses the original loaded full inertia and incoming momentum mismatch,
prepares the native prescribed-base row mobility, and solves all row reactions
together. One final native impulse event applies the capture mismatch and
`A^T lambda` simultaneously. The actual endpoint retains the original speed
and representability admission; an inadmissible unconstrained intermediate
cannot reject a valid constrained result. Passive impulses do not consume
motor force/torque ratings.

The result wraps the original capture with row impulses, multipliers, initial
and final row residuals, conditioning, and independently bounded midpoint
reaction work. The native support impulse and signed support work remain
external. Actual world linear/angular momentum and RNEA-projected joint momentum
are checked against the support and declared ideal reactions respectively.
The same private lifecycle serves the original unconstrained API, and empty
rows reproduce that capture exactly.

Capture retains its nonnegative complete relative-motion quadratic Q. Its
identity is `delta_K = W_support + J_constraint dot v_final - Q`.
For exactly tangent motion the last pairing vanishes. For finite admitted row
defects, the final-rate pairing remains visible in the existing energy residual,
while actual midpoint constraint work has its own cap. Neither term is silently
subtracted from Q, double-counted as loss, or assigned to a thermal body.

Seven focused tests use independent COM/spin and reduced-inertia calculations
for geared finite recoil, positive/negative moving-support work, a scalar-root
lock, full-tensor branched motion, frame and inertial-boost covariance, finite
tangency work, row scaling, actual endpoint limits, budget refusal and cancellation.
The continuation case captures a payload, advances the loaded geared mechanism
through the actual supported implicit path, and releases its exact inventory.
This admits existing constant affine gears/locks, not constraint activation,
backlash, varying ratios or nonlinear point loops.

## Supported constrained payload release and ejection

`articulated::motion::inventory::supported::constraints` now exposes passive
`release_payload` and driven `release_payload_with_impulse` alongside capture.
They share the original inventory identity, rigid payload state reconstruction,
physical point/couple convention and common release lifecycle. Active release
uses the exact restored retained inertia for the native ideal-row mobility,
then applies the separation impulse and `A^T lambda` in one final native
prescribed-base event. The outgoing body retains its original full tensor and
inverse-inertia response.

Actual world momentum changes are audited against the independent support
impulse; native joint momentum changes also include the declared ideal row
impulse. The energy identity is
`delta_K = W_support + W_separation + W_constraint + defect`.
The separation pair begins on one rigid velocity field, so its initial work
cancels and its source work comes from the actual velocity jumps. Constraint
midpoint work and signed support work remain separate. The final balance
retains compensated reaction products before rounding public work accounts.
Source and absolute reaction-work budgets are independently enforced; the
existing source-only policy refuses negative separation work. No capture
projection loss or thermal sink is introduced into release.

Passive and zero-impulse release admit initial row residuals and the retained
inertia's normalized mobility without projecting rates or applying a reaction.
They retain exact original state and physical balance. Numerical probe rates
do not receive physical speed limits; the completed native event does.

Five further test groups extend the shared constrained inventory target to
12 cases. Independent retained-inertia reductions check geared recoil and
signed support work; COM/spin and a separate centroidal solve check full-tensor
payload response and off-centre pair work. Other cases retain finite reaction
work, source/rate/rank refusal, passive exactness, frame/boost covariance and
cancellation. The continuation case now covers capture, actual constrained
implicit motion, and ejection with the advanced state and original attachment.
Existing constant affine rows are supported; nonlinear loops, changing ratios,
contact-triggered release and finite-duration actuator models remain separate.

## Simultaneous point ejection at prescribed relative velocity

`assembly::release::eject_point_masses` closes the finite point-release model
with prescribed FINAL parcel velocities relative to the recoiling retained
rigid body at each physical attachment. All jets use one original body snapshot.
For original total mass T, full central inertia I, arms r_i and relative
velocities c_i, conservation gives dv=-sum(m_i c_i)/T and
dw=-I^-1 sum(r_i cross m_i c_i), hence J_i=m_i(c_i+dv+dw cross r_i).
The native `release_point_masses` owner then admits the actual retained
COM/inertia, applies these impulses and computes their complete joint work.

The projection uses O(inputs) fixed-order compensated reductions and one native
original-inertia mobility action, without a dense 3N constraint solve. Its
translation chart is anchored at the largest parcel's requested velocity and
uses the represented retained-mass fraction directly; near-total central
ejection therefore retains the small reduced-mass impulse. Computed native
impulse descriptors reuse the staging allocation in the returned result.
Every descriptor consumes its part slot, and every nonzero SOLVED impulse
consumes an impulse slot. Even a zero-relative-speed parcel can need impulse
to follow the shared recoil induced by another jet.

Explicit finite nonnegative allowances bound supplied separation work and the
maximum component of actual relative-velocity error. The original momentum and
energy gates remain active. The velocity residual uses actual canonical
momenta in one common translating chart and actual returned parcel/body
positions; rounded-away recoil or geometry cannot be hidden by intended
velocities. The source owner must debit the accepted supplied work. Cancellation
or any late failure publishes no changed body, parcel or source state.

Six focused G0/G3 regressions specify central reduced-mass response and extreme
retained-mass ratio, parallel/opposed coupled jets, a full-tensor asymmetric
case whose zero-speed parcel still needs impulse, covariance and subdivision,
large common boosts and actual-velocity refusal, passive native parity, source
work and input budgets, and late cancellation. Pinned Rustfmt parsing and
independent rational-oracle checks passed; native compilation and execution
remain unrun in this session. This is instantaneous point detachment, not a
continuous nozzle/rocket law, pressure-thrust model, finite-fragment release,
fracture model or propellant caloric source.

## Atomic finite point-material transfer between moving bodies

`assembly::transfer::transfer_point_masses` connects the native simultaneous
release and capture owners in one atomic local event. Donor descriptors and
receiver attachment arms use their respective original COMs and body axes.
Both bodies receive full central-inertia updates and shared recoil. No free
parcel endpoint is published alongside the receiver, preventing duplicate
material inventory at this interface.

Attachment admission checks origin-relative geometry before forming large
absolute world positions, and also checks the represented release/capture
locations. This is instantaneous handoff at a common physical point, not
teleportation across a gap or a free-flight/contact-discovery law. The caller
still owns geometry membership, material identity and the choice of event.

`PointMassMomentumCapture` and `capture_point_masses_with_momenta` add direct
canonical momentum ingress to the same batch projection used by authored
velocities. They do not divide momentum by mass and multiply back to form the
receiver's momentum. Fused subtraction of the common boost supplies relative
momenta for spin and loss. The original velocity API and singleton fast path
remain available; zero canonical mass requires zero momentum.

A combined gate reconstructs actual two-body mass, linear momentum, spin plus
orbital angular momentum about the declared original donor COM, and kinetic
change minus separation work plus sticking loss. Polarized changes and a
common translating frame retain small exchanges under a large uniform boost;
actual first-moment changes remain in the angular balance. Separation work
requires a caller-owned source, and sticking loss a caller-owned destination.
Neither is silently deposited as heat or given thermodynamic authority.

The work cap counts both original bodies and every descriptor, including
zeros. Positive capture impulses and nonzero separation impulses share one
aggregate impulse cap. The O(N) native passes and staging remain bounded;
each phase admits one aggregate inertia, with no per-parcel factorization.
All refusals and cancellation leave both borrowed bodies unchanged.

Seven focused regressions cover canonical handoff, independently known
finite-body recoil/work/loss, full tensors and independent rotating frames,
common boosts, batch permutation/subdivision, attachment gaps hidden by large
coordinates, combined angular rejection, zero identities, invalid inventory,
global budgets and late cancellation. Native execution is required before
claiming numerical conformance; formatter parsing is not that evidence.

## Balanced physical nozzle arrays and full spatial exhaust transport

`rocket::rigid::nozzles` prepares simultaneous central nozzles and opposite-arm,
equal-rate point-nozzle pairs. Relative exhaust velocities belong to the LOCAL
rotating port. Structural pairing makes the outgoing mass first moment zero,
matching the existing fixed-COM model without a fitted balance tolerance.
From the same rates/arms/velocities it derives the resultant reaction force,
relative reaction couple, full co-rotating exit inertia, and local-relative
kinetic power. Zero/rank-deficient exit moments are physical and admitted by
positive point-mass construction; no diagonal approximation or artificial
positive inertia is inserted.

`PreparedNozzleArray::forcing` maps these physical moments into the original
full-tensor forcing. The thrust arm is zero and all relative tangential jet
momentum is exported through the original swirl channel, so pure couples are
supported and relative torque is counted once. `TensorRocketState::advance_nozzles`
reuses native mass/inertia, analytic translation, bounded rotational substeps and
exact burnout segmentation. Additional loads use the same prepared forcing with
the existing `advance_spatial`; neither propulsion nor its moment is duplicated
as an external wrench. Rates/geometry are held during a call and all valves close
at the original propellant floor.

`flux_at` reports instantaneous physical mass, absolute world linear momentum,
angular momentum about the COM, and angular momentum about a supplied world
point. It includes COM velocity advection, the full co-rotating orbital spin of
the exported particles, relative jet angular momentum, and the COM orbital term.
It is an instantaneous transport query, not an independently integrated impulse
receipt. Local-relative kinetic power is a diagnostic only: no extra source-work,
heat, pressure-work or caloric debit is introduced.

Preparation is allocation-free, bounded to 16,384 groups, compensated in supplied
order, and cancellation-polled in both passes. Stepping/query storage and work
are independent of nozzle count. Focused G0/G1/G3/G4 tests cover physical
individual-particle flux sums, origin/Galilean covariance, pure-couple depletion
and burnout refinement, rotating off-axis venting, the original analytic central
rocket, and admission/cancellation. Arbitrary unbalanced outflow, moving COM,
internal relative-flow momentum, changing reservoir geometry, pressure thrust,
thermal/chemical plume energy and complete open-system passivity remain outside
this bounded model. The reservoir inertia law remains an independent physical
declaration; balanced nozzle geometry does not determine it.

## Simultaneous closed finite-body material-transfer networks

`assembly::transfer::network::transfer_point_mass_network` composes the original
native batch release and canonical-momentum capture owners across a directed
graph. Every outgoing parcel uses its donor's ORIGINAL snapshot. All releases
complete before any capture, so cycles and nodes that both give and receive
never re-export newly received material during the same event. Each receiver
captures its whole incoming batch together after rebasing original attachment
arms by its own release COM shift.

Two prefix-indexed adjacency tables route parcels in O(bodies + edges) work
and storage, retaining original body and per-body input order. Reused descriptor
buffers avoid per-edge staging allocations. Each participating body admits at
most one release and one capture inertia; no dense network mobility matrix or
repeated pairwise recoil solve is introduced. Canonical parcel momenta are
consumed exactly once into the receivers, and no free parcel inventory is
published alongside the final bodies.

Original local attachment geometry and represented post-release positions
must both coincide. Each intermediate retained donor must keep positive mass
and admissible inertia even if incoming material later replaces its inventory.
The event is instantaneous physical handoff, not finite flight, contact discovery,
material/species mixing, or a continuous-flow integration scheme.

The original shared controls count all bodies plus every edge, including zero
mass. Positive capture and nonzero separation impulses share the aggregate
impulse cap. Native per-phase balances and a common-origin global mass,
linear/angular momentum and polarized kinetic/work/loss gate admit only the
complete event. Unchanged isolated bodies contribute an exact zero change,
without evaluating unnecessary absolute energies or distant moment arms.
Separation work and sticking loss remain separately owned per body and globally;
source debits and thermal destinations are explicit caller responsibilities.
Cancellation and late failure publish no endpoint and preserve all inputs.

Six focused regressions cover an independent three-body circulation oracle,
shared incoming recoil, a rotating full-inertia off-center cycle with known
central inertia, edge/body permutations and parcel subdivision, pair parity
with canonical momentum and isolated owners, large common boosts, budgets,
invalid intermediate inventory and terminal cancellation. Existing pairwise
transfer reuses the same generalized physical-balance and attachment helpers.
Native compilation and numerical execution remain necessary for conformance;
formatter parsing and analytical checks do not substitute for them.

## Coupled finite-body revolute hinges

`rigid_rattle::distance::graph::advance_articulated_graph` extends the existing
native coupled graph with `RevoluteJoint`, alongside its original spherical
and positive-distance joints. A hinge contributes three coincident-point rows
and two body-relative alignment rows, leaving relative axial rotation free.
The two axes are normalized once; a transverse basis is prepared in body two.
The positive directed-axis branch is explicit, and opposite axes refuse.

For unit world axes `a` and `b` and body-two transverse directions `t_i`, the
angular equations are `g_i = L * a.dot(t_i)`. The exact torque Jacobian on body
two is `L * (t_i cross a)`; body one receives its negative. Thus both axes follow
their evolving native attitudes, and the physical reaction is a pure couple.
The complete force/couple `J M^-1 J^T` includes all cross responses between
coincident-point rows, angular rows and shared full-inertia bodies. It is the
original bounded dense owner, not sequential independent joint solves or a
redundant six-row hinge.

The caller supplies a positive `angular_length_m = L`. Original graph position
and velocity tolerances `(p,v)` remain in metres and metres per second. The
angular admission is `sin(axis_angle) <= p/L` on `a.dot(b)>0`, with the Euclidean
norm of the two alignment derivatives at most `v/L` per second. At exact
alignment that norm is the relative tilt speed; away from exact alignment it
is the derivative of the stated equations, including their moving basis.
`L` must exceed `p`. Changing `L` without changing the controls changes these
authored angular tolerances. The complete hinge couple at each impulse stage
is bounded by `L * maximum_impulse_ns`; each coincident-point reaction retains
its separate original vector impulse cap. Bounds use the actual summed couple,
including finite-residual nonorthogonality, not independent component caps.

Every pure couple calls the original `apply_momentum_impulse` on both bodies.
`JointReactionImpulse` distinguishes untouched native point-pair receipts from
`CoupleImpulseReceipt`, whose two native updates retain requested impulses,
actual representable increments and residuals. Paired angular work uses the
commanded world couple and double-double midpoint relative angular velocity.
Paired kinetic change uses the actual body-momentum increments. Their difference
and actual total angular increment remain visible; no fictitious force arm or
pose repair hides them. The joint receipt also retains original free steps,
point-gap angular impulse, scalar residuals and authored-joint row spans.

The shared rigid impulse owner evaluates fixed-mass/principal-inertia storage
from compensated endpoint momenta, including non-dyadic mass and inertia
quotients. Native endpoint-velocity work stays independent: both represented
impulse loss and the velocity-rounding defect remain in its work residual.
Point pairs reduce the original storage and commanded point-work product layers
before rounding, rather than adding rounded per-body energies or residuals.
`RigidKineticEnergyChange` exposes this allocation-free storage reduction to
coupled consumers; each body must enter once with fixed mass and material-frame
inertia, and the caller retains its own work, heat, and cancellation ownership.

The generic `DistanceGraphStep` defaults to its original point receipt, so
existing distance, driven-distance and spherical APIs retain their event type
and native event layout. The current geometry-cache cost is specified in the
moving-plane/prismatic extension below. Shared trial buffers and mobility
factorization
serve both paths. An articulated run adds one fixed 32-row stack bank and one
bounded joint-row map; no per-trial allocation or per-step axis preparation is
introduced. A hinge consumes five of the existing 32 scalar rows. Native body
travel, local solve/backtrack/rank and whole-transaction memory/work/cancellation
bounds remain in force. Any refusal publishes no partial state.

Nine focused G0/G1/G3/G4 tests cover analytical unequal-inertia hinge-chain
torque and mass recoil, free axial spin, an independent continuous two-rotor
gyroscopic precession reference with refinement/reversal, long off-center
full-tensor mixed-joint momentum/work balance, frame/boost/endpoint/order and
explicit-scale covariance, point-only native parity, opposed/invalid axes and
large-boost admission, complete vector caps, the actual connected 32-row
boundary, rank/payload/trial limits and every observed cancellation checkpoint.
All 24 pre-existing distance, driven and spherical graph tests also pass in the
bounded actual-source native build. This evidence does not substitute for
unrun full-workspace/DSR gates or independent experimental validation.

This remains a local smooth bilateral hinge between finite rigid bodies.
Joint friction, motors, angular stops, backlash, impacts, fixed supports,
arbitrary dependent constraints, nonlinear state enclosures, exact energy,
global branches, symplecticity and the other integrator lanes remain unclaimed.

## Finite dry friction on a native revolute joint

`rigid_rattle::distance::graph::revolute_friction::apply_revolute_friction_impulse`
applies an axial drive/friction event to one existing hinge while every original
bilateral graph row reacts. The complete full-inertia constrained response is
`a M^-1 a^T - b^T (J M^-1 J^T)^-1 b`, with `b = J M^-1 a^T`.
The existing geometry, mobility factorization and native point/couple impulse
owners preserve shared-body recoil, including off-centre hinges and connected
spherical, distance or prismatic joints. A locked or unresolved axial response
refuses without adding compliance.

The caller declares physical torque-impulse capacities. An already moving
hinge uses only kinetic capacity, including when it arrests. A separate
`RestingDrive` phase requires exactly zero initial axial relative rate and
applies the drive and friction together: a subthreshold drive is held with
exactly zero resultant impulse and work; breakaway uses kinetic capacity.
No numerical speed threshold spends static capacity on moving spin.

Actual native states retain their original poses. The report separates signed
drive work, friction work, bounded bilateral residual work and independently
reconstructed kinetic change. Native total linear/angular momentum increments
and the actual outgoing axial rate must meet their authored gates; neither
momentum nor small arrest residuals are projected away. Cancellation, effort,
memory, rank or conservation refusal publishes no candidate.

This is one frictional joint at a fixed-pose impulse event, with at most the
original 32 bilateral rows. It reuses the bounded dense graph factorization
and shared Coulomb impulse projection. The caller owns capacities, phase
selection and subsequent native time integration; simultaneous friction at
several joints, continuous stick-slip transitions, bearing load/radius models,
elastic preslip and heat allocation are outside this operation.
`tests/rigid_revolute_friction.rs` covers analytic rotors, static hold and
breakaway, off-centre and third-body recoil, full-tensor frame covariance,
native RATTLE continuation and atomic refusal.

## Coupled finite-joint friction and native slider phases

`rigid_rattle::distance::graph::joint_friction::apply_joint_friction_impulses`
solves explicitly supplied axial friction capacities on several original
`ArticulatedJoint::Revolute` and `ArticulatedJoint::Prismatic` joints at once.
It retains the full projected port mobility
`P M^-1 P^T - B^T (J M^-1 J^T)^-1 B`, including every cross response through
the existing bilateral graph. One original bilateral factorization supplies
all port reaction columns. Hinge coordinates use the original authored length
`L`: rate `L * relative_spin`, impulse `angular_impulse / L`.
The rail coordinate reuses the original moving-plane Jacobian along the
carrier's joint-frame Z axis, including its rotating-axis/current-point terms
and reciprocal force arms. Full body inertia and finite carrier recoil remain
in the native response.

The shared `fs-contact` bounded friction owner solves the complete port problem.
Native point impulses, couples and summed bilateral reactions then act on the
original bodies in increasing joint-index order. Actual final port rates,
original bilateral closure, capacity support, native kinetic change, and
linear/angular momentum increments are independently gated before returning the
candidate. Signed midpoint work at one port may be positive through coupling;
only the aggregate friction work is a mechanical-loss statement. Bilateral
residual work is separately bounded and neither quantity allocates local heat.

`apply_prismatic_friction_impulse` additionally provides explicit single-slider
`Sliding` and `RestingDrive` phases. An exactly resting drive below its static
capacity has exactly zero resultant impulse and work. Breakaway uses the
kinetic capacity; pre-existing slip cannot borrow the static capacity or be
captured by a numerical speed band. The existing single-revolute API retains
its physical behavior and shares the same native momentum-accounting helper.

These are fixed-pose events with caller-supplied impulse capacities, not a
continuous trajectory, inferred load law or automatic multi-joint friction
history. Up to six selected hinges/sliders fit the original 32 bilateral rows;
other graph constraints share that same row budget. Port and cross-response
matrices use bounded stack banks. Empty port sets preserve the admitted native
state without factoring a friction mobility. Duplicate, locked/dependent,
numerically non-SPD, nonrepresentable, cancelled, over-budget and unresolved
events refuse without changing borrowed states. Numerical effort bounds never
replace or shrink the authored physical friction capacities.

## Finite-time friction in native articulated mechanisms

`joint_friction::relaxation::advance_joint_friction` evolves original hinge and
slider velocities over a declared duration under held axial torques/forces and
Coulomb capacities. It consumes the shared continuous `fs-contact` flow with
equal static and kinetic bounds, including coupled slowing, isolated arrest,
holding and renewed motion. Successive declared load segments may change the
drives or capacities. The entire original graph participates through one
bilateral factorization and the full projected port mobility; the exact moving
rail Jacobian, original full inertia and finite carrier recoil are preserved.

Preparation and native impulse application are shared with the instantaneous
joint-friction owner. The physical acceptance tests remain distinct: continuous
flow checks its whole velocity path and work, while the impulse API retains its
endpoint support/KKT law. Integrated force/couple actions update the original
momenta, then actual native endpoints are checked against both the original
mobility response and the resolved flow endpoint. Original bilateral closure,
linear/angular momentum, native kinetic change, physical capacity residuals and
continuous drive/friction work must all pass before a candidate is returned.

Signed and absolute hinge rotation or slider travel are physical path
integrals. Each continuous friction port has nonpositive path work, even when
its net impulse and the two outer endpoint velocities would suggest otherwise.
Drive work is integrated over that same path; bilateral residual work remains
separate. Floating capacity/drive-conversion residuals are reported and bounded,
never clipped into authored actions or relabeled heat.

Geometry, inertia and constraint Jacobians are fixed during this phase, and
poses are unchanged. The surrounding trajectory owner must transport poses and
re-admit changed geometry. Normal bearing loads, different static/kinetic
history, wear and thermal allocation are not inferred. The existing six-port
and 32-bilateral-row limits, native resource gates and atomic cancellation
semantics remain in force. Existing single-hinge, single-slider and simultaneous
impulse APIs retain their original physical laws.

## Driven distance actuators inside finite ball-and-hinge mechanisms

`rigid_rattle::distance::graph::advance_driven_joint_graph` accepts a single
`DrivenJoint` sequence containing quadratic distance actuators, spherical
joints and revolute hinges. Every scalar row enters the original simultaneous
full-inertia force/couple mobility and native impulse/free-step/impulse solve.
An actuator therefore recoils every connected finite body and transmits the
actual off-center forces and hinge couples. No part of the mechanism is pinned
or advanced by a separate integrator or sequential joint solve.

The original `DistanceDrive` command owner prepares each radial row's positive
length, initial/final range rate and representation residuals. It still checks
interior collapse and unrepresentable nonzero displacement. Cartesian ball
and angular hinge rows retain their original geometry and zero prescribed
rate. Default distance commands supply fixed rods. The rate owner cancels
common COM translation before adding spin, including for default commands;
rounded laboratory point velocities cannot hide a spinning attachment.

Every original native point/couple event remains in the result. Actuator
source work is its two scalar impulses paired with the corresponding commanded
range rates. Passive joints contribute exactly zero source work. Their actual
ordered force/couple work, including shared-body cross terms, remains in the
constraint-work defect. The result also pairs total reaction work independently
and reports the requested residual angular impulse of finite point closures.
Native momentum/kinetic-work representation defects remain in the underlying
receipts; no energy discrepancy is erased or labeled as physical heat.

One continuation vector carries each updated command, per-joint work,
representation residual and scalar-row span. The original 32-row bound counts
one row per rod, three per ball and five per hinge. Two fixed stack banks extend
the original reused native trial buffers; no per-trial allocation is added.
Live-payload admission includes the complete continuation vector. Rank, local
branch, whole-vector impulse, angular travel, nonlinear trial and cancellation
limits still refuse atomically. The shared internal free-motion seam uses the
actual native free-step receipt type and includes caller-owned load-bound
storage in the same budget; the public held-wrench path retains its native map.

Seven focused G0/G1/G3/G4 tests cover analytical finite group recoil and hinge
torque, a free-floating actuated hinge with a finite spherical payload against
an independent conserved-angular-momentum quadrature and refinement, long
off-center full-tensor work/momentum balance and reversal, native passive and
rod-only parity, frame/boost/endpoint/order covariance and hidden-spin admission,
the actual connected 32-row boundary, command/rank/trial/payload refusals and
every observed cancellation checkpoint. All 59 prior native MBD regressions
also pass in the actual-source bounded native build. Full-workspace/DSR gates
and independent experimental validation are not substituted by this evidence.

This remains an ideal bilateral distance actuator with an external energy
source. Motor efficiency/storage, commanded hinge spin, friction, stops,
impacts, dependent-row mechanisms, global branches, exact energy and
symplecticity remain outside the claim.

### Spatial force moments through complete finite-body mechanisms

`advance_spatial_articulated_graph` and `advance_spatial_driven_joint_graph`
accept the original `PreparedSpatialLoads` for every finite body in the passive
or actuated rod/ball/hinge graph. Each
nonlinear trial calls the original `RigidBodyIntegrator::step_spatial`: a
world-held point force follows its rotating material arm, while body-held and
world-held couples retain their declared frames. The complete native
`SpatialStepReceipt` is returned for every free step, including source impulse,
actual work and its finite-step defect. Joint reactions retain their original
finite recoil, point-gap angular impulse and force/couple work accounts.

Pair and graph paths share one statically selected native free-motion adapter.
The force/couple mobility, native full-inertia response, two impulse stages,
nonlinear trials, rank checks and authored-joint mapping remain the same
owners. `DistanceGraphStep` and the joint receipts default to their original
free-step type, preserving the held-wrench API and native free-step layout. The
spatial graph budgets one torque-bound scalar per body in addition to its actual
larger native receipts and joint map. Each conservative torque bound is prepared
once from the fixed-size load representation; no load-set traversal, fresh
inertia factorization or allocation occurs within a trial.

The spatial driven entry uses the same target bank, quadratic command checks,
per-joint continuation and source-work owner as the held-wrench driven entry.
External spatial work, prescribed actuator work and passive constraint defects
remain separate. Passive-only motion has identical native receipts; default
radial commands use the existing compensated rate arithmetic and agree with
the legacy static-rod result at roundoff rather than promising bitwise identity.

Seven groups in `tests/rigid_spatial_graph.rs` check an independently integrated
nonlinear loaded hinge chain with second-order refinement and reversal,
off-center full-tensor mixed-joint source/constraint momentum and work closure,
world-frame and common-boost covariance, exact native free-body/source-receipt
identity, a loaded actuator with independent unequal-mass recoil, rotor dynamics
and signed actuator work, complete payload admission and every observed
cancellation point. All 73 combined tests (the 59 earlier pair, graph,
driven-distance, ball, hinge and native spatial-load tests plus fourteen new
mixed-drive and spatial groups) pass in the bounded actual-source build. This is not a full-workspace
or DSR gate. Held loads, resolved local travel and the original smooth bilateral
joint assumptions remain required; arbitrary state-dependent force laws,
contact, global branches and exact energy conservation remain outside this map.

## Coupled moving planes and finite prismatic mechanisms

`PointPlaneJoint` and `PrismaticJoint` enter the same native graph through
`ArticulatedJoint` and `DrivenJoint`, including both original spatial-load
entry points. They can share arbitrary finite bodies with rods, balls, hinges
and quadratic radial actuators. Every point-force and angular-couple cross
response enters the existing complete full-inertia mobility and simultaneous
bounded impulse/free-step/impulse solve. No carrier is pinned or treated as
infinitely massive, and no pose correction or sequential subsystem solve is
introduced.

`PointPlaneJoint` reuses the original finite point/plane geometry owner. Its
plane normal is normalized once and only the body-normal anchor offset is
retained. The instantaneous equation uses the material point relative to the
carrier COM and the carrier's current world normal. Its rate includes the
carrier angular velocity at that CURRENT sliding location, not at a fixed
authored anchor. Native equal/opposite impulses act at the material point and
its projected point on the plane. The cached body-local carrier arm therefore
changes as the point slides; mobility, trials, endpoint reactions and receipts
all use that same actual arm. Both original pair and new graph paths share
this preparation and geometry implementation, including compensated relative
positions and COM/spin rates, irrelevant tangential-anchor removal and reaction
arm representability admission.

`PrismaticJoint` contributes two such transverse plane rows followed by three
relative-orientation rows. Local joint Z is the free sliding direction. The
two caller-supplied normalized frames map joint axes into the respective native
COM body frames, so neither inertia principal axes nor matching body attitudes
are assumed. The complete joint frames remain aligned while arbitrary axial
translation is permitted. The two planes are prepared once from the carrier
frame and share the original moving-plane torque and finite recoil treatment.

For current world joint frames `Q0` and `Q1`, let the relative quaternion be
`Q0^-1 Q1 = (w,v)`. The angular gap vector is `g = 2 L w v`, where the caller
supplies positive `angular_length_m = L`. Its exact world Jacobian for component
`i` is `n_i = L Q0 [w*w e_i - v_i v + w (v cross e_i)]`; the rate is
`n_i.dot(omega1 - omega0)`. These analytic moving-frame Jacobians apply native
equal/opposite pure couples through the original momentum-impulse owner.
Production uses no angular finite differences or additional inertia
factorization. Common angular motion is cancelled before compensated pairing,
and all reaction work and representability residuals remain in native receipts.

The admitted orientation chart requires `w*w > 1/2`, equivalently relative
rotation magnitude below pi/2. This excludes the sine chart's half-turn false
zero and quarter-turn rank loss. `L` must exceed the absolute position tolerance
`p`. The complete transverse gap/rate norms obey the original `(p,v)` controls;
the complete angular gap norm is `L*abs(sin(theta))`, and its derivative norm
obeys the same speed tolerance. At alignment the latter is
`L*norm(omega1 - omega0)`. The complete two-row transverse impulse is bounded
by `maximum_impulse_ns`. The summed physical three-row couple has its separate
bound `L*maximum_impulse_ns`, including finite-residual nonorthogonality. Scalar
component caps cannot replace these vector bounds. Native angular travel,
moving-normal travel, rank, backtracking and endpoint admission remain binding.

A rod or standalone point/plane consumes one scalar row, a ball three, and a
hinge or slider five of the existing 32-row bound. The original pair and graph
now share one `Geometry` representation containing both actual reaction arms.
The three graph geometry banks therefore add `3*2*sizeof(Vec3)` logical bytes
per scalar row (144 bytes for the native three-f64 vector), including legacy
rod/ball entry points. Exact current native receipt, row-map and continuation
sizes enter the existing preallocation byte admission. The bounded row banks
remain on the stack, and the existing trial buffers are reused without
per-trial allocation. Original distance/spherical native event types remain
unchanged. Every failure and cancellation leaves all borrowed inputs intact.

Three groups in `tests/rigid_plane_graph.rs` cover an independent coupled
shared-carrier acceleration solution with original pair parity; two sliding
masses exchanging orbital momentum with a rotating finite carrier against an
independent continuous ODE, second-order refinement and reversal; and off-center
full-tensor moving planes combined with radial actuation and spatial forces,
with momentum/work closure, covariance, permutation, limits and cancellation.
Six groups in `tests/rigid_prismatic_graph.rs` cover exact finite rail recoil,
transverse load sharing, axial torque and free travel; rotating-rail centrifugal
motion against an independent continuous ODE with refinement/reversal; a driven
closed slider-crank against independent finite-carrier angular-momentum
quadrature and actuator-work convergence; a full-tensor rail/hinge/actuator with
native spatial receipt identity, repeated linear/angular momentum and work
closure, rotation/common-boost covariance; a genuinely connected 32-row graph,
exact byte boundary, frame-branch/hidden-rate/vector-cap refusals; and every
observed cancellation checkpoint with identical retry. A native inline unit
test independently derives orientation gaps from transformed axes and checks
all moving-frame Jacobian components against central rotational differences.

These nine new integration groups and the independent Jacobian unit test pass
in the bounded actual-source native build, alongside all 73 previous affected
native regressions. This is numerical source-level evidence, not a full
Cargo workspace or DSR gate or independent experimental validation. The wider
multibody work remains open. The extension admits smooth bilateral finite rails
and infinite planes, with explicit local orientation and independent-row
requirements. Travel stops, friction, unilateral impact/contact, rail motors,
arbitrary dependent constraints, global branches and exact energy/symplectic
claims require their own models.

## Finite revolute angular actuation in coupled native mechanisms

`RevoluteMotorJoint`, `AngularDrive` and `DrivenJoint::RevoluteMotor` prescribe
the remaining axial angle of a finite hinge in both `advance_driven_joint_graph`
and `advance_spatial_driven_joint_graph`. Two material points coincide, and the
current first joint frame turns about its local Z by
`delta(t) = rate*t + acceleration*t*t/2`. Body two follows that moving frame.
Both bodies retain finite mass, full intrinsic inertia and recoil; the joint
frames need not coincide with either inertia's principal axes. Motors share
the same coupled solve with other motors, radial drives, rods, balls, passive
hinges, moving planes and sliders. An ideal zero-rate motor holds its relative
angle and has six rows; the passive five-row hinge still rotates freely.

Each motor contributes the original three Cartesian point closures and three
full-frame angular rows from the prismatic geometry owner. For relative joint
quaternion `(w,v)`, `g = 2 L w v` and the existing exact world Jacobian is
`n_i = L Q0 [w*w e_i - v_i v + w (v cross e_i)]`. With current world motor axis
`a = R_body0 frame0 Z` and commanded relative spin `nu`, its hidden rate is
`n_i.dot(omega1 - omega0 - a*nu)`. The command is subtracted in compensated
arithmetic before pairing with the Jacobian. The advanced first frame enters
endpoint geometry directly; no finite-difference angular Jacobian, new free
integrator, pose repair, fixed support or sequential pair solve is introduced.
The existing full-inertia force/couple mobility retains all shared-body cross
responses, nonlinear backtracks and actual native endpoint impulses.

The quadratic angular increment and final rate use the same compensated
polynomial owner as radial actuation. Admission checks the maximum absolute
command excursion over the complete interval, including its interior turning
point, against twice the native per-body angular-travel cap. A hidden full turn
or excessive out-and-back command therefore refuses even if its endpoint frame
would coincide with the start. Actual native body travel, the positive-cosine
orientation chart, scaled pivots and complete vector impulse bounds still apply.
Moving-normal admission includes the separately bounded commanded-frame travel.

The original native quaternion exponential prepares the advanced local frame.
Compensated relative-quaternion products reconstruct its axial phase and any
unwanted transverse tilt, including canonical sign changes across accumulated
half turns. `AngularDriveStep` exposes the represented displacement, axial
representation residual, transverse frame residual, continuation-rate residual
and complete interval excursion. `L*(abs(phase residual) + transverse tilt)` must
fit the position tolerance, as must its sum with the final angular gap norm.
The represented continuation rate and actual advanced frame are also checked
against the hidden-rate tolerance. An unrepresentable command or nonzero
displacement that vanishes in its frame refuses. Continuation carries the
advanced frame and final rate, allowing repeated revolutions without recovering
an absolute wrapped angle from a body pose.

Signed motor source work pairs each ACTUAL native world couple impulse with
the prescribed angular velocity of the current motor axis at its impulse
stage: `couple_on_second.dot(a)*nu`. The source account therefore includes full
tensor and shared-body effects from the coupled solution. Positive values
supply energy; negative values recover it. Original native paired reaction work
minus this source remains `constraint_work_defect_j`, including all ordered
cross terms and passive point/angular work. Passive joints supply zero actuator
work. Prepared spatial-load source work and native free-step numerical defects
remain separate, and all original momentum and impulse receipts are retained.
No defect is interpreted as physical heat or discarded to enforce energy.

The motor's six actual scalar rows count toward the existing 32-row limit.
Its three point impulses and summed three-row physical couple have the same
separate complete-vector caps as the original point and full-frame owners; the
couple cap is `L*maximum_impulse_ns`. Existing fixed row/target banks stay on the
stack, current native continuation sizes enter live-byte admission, and original
trial buffers are reused with no per-trial allocation. Every refusal and
cancellation publishes no partial state or continuation and leaves borrowed
inputs intact. Angular command failures identify the first motor angular row.

Six integration groups in `tests/rigid_angular_drive_graph.rs` cover a loaded
three-rotor chain against independent analytic angular-momentum and torque
balances, signed supply/recovery, nine-radian continuation and reversal; two
full non-diagonal tensors against an independent conserved-angular-momentum
quaternion ODE with second-order state and source-work refinement and reversal;
an off-center motor/slider/radial mechanism under native spatial loads with
actual couple-work pairing, free-source receipt identity, repeated P/L/work
closure and rotation/common-boost covariance; a connected 32-row motor chain,
exact live-byte boundary and complete physical couple cap; full-interval winding,
turning-point and hidden-rate admission; and every observed cancellation point
with identical retry. All six new groups and 82 previous affected integration
tests pass in the bounded actual-source native build. This is numerical
source-level evidence, not a full Cargo workspace/DSR gate or experimental
validation. These are prescribed ideal bilateral actuators; finite motor
power/storage, efficiency, friction, stops, unilateral contact, arbitrary
dependent constraints and global branches still require their own models.

## Finite rigid separation at prescribed relative linear and spin velocities

The native assembly owner exposes
ejection::RelativeSeparationTwist and RigidAssembly::release_relative_twists.
Each component has a prescribed final COM velocity relative to a chosen
reference component's FINAL rigid velocity field at that COM, plus a prescribed
angular-velocity difference. Both vectors use the original assembly axes.
The reference entry is exactly zero. All full intrinsic tensors, distinct
component axes, COM arms, and the original composite inertia participate.

One shared recoil solves linear and spin-plus-orbital momentum together:
dv = -sum(m_i c_i)/M and
dw = -I^-1 sum(r_i cross m_i c_i + J_i eta_i).
The implementation subtracts an equivalent common rigid twist before fixed-order
compensated reduction. The largest-mass component anchors linear data and the
largest intrinsic principal moment anchors angular data. Small reactions can
therefore survive an inventory or spin inertia dominated by one component.
This is O(parts) work/storage with one cached native angular-mobility action,
not a dense 6N constraint solve or independent sequential pair impulses.

Equivalent star force/couple pairs against the reference enter the ORIGINAL
native release owner. Point-force moments and intrinsic spin couples are
distinct. The native aggregate per-body impulse work includes all cross terms,
and every released body retains its original material coordinates and full
mass properties. Zero commanded relative motion can still need impulse to
follow the common recoil. Only nonzero solved pairs consume the impulse cap.

The returned canonical momenta and represented COM poses determine the actual
relative linear and angular rates. A common translational chart avoids hiding
small velocity changes behind separately rounded laboratory velocities.
Native momentum/work admission, both relative-rate tolerances, the available
source-work limit, and final cancellation must all pass before publication.
The source must explicitly debit the accepted mechanical work. Contact,
actuator placement, fracture, finite-time forcing and thermal ownership remain
caller responsibilities; the star wrench decomposition does not infer them.

AssemblyPartition::release_relative_twists applies the same law to cached
retained groups. It uses group indices and full group inertia, while the part
cap still counts all original components. It never releases leaves and then
re-latches them, avoiding fictitious sticking loss and incorrect separation work.

Six focused G0/G3/G4 regressions cover independent linear/angular reduced-inertia
and spin-orbit oracles, full-tensor frames and reference/component permutations,
retained-group rigidity, extreme inventory/inertia ratios, source/budget/input
refusal, unrepresented actual endpoint rates, and terminal cancellation.
Pinned formatter parsing and source/algebra checks are not native numerical
conformance: compilation and these tests remained unrun in this environment.


## Mixed spherical scenes with accelerating plane boundaries

`articulated_contacts::trajectory::mixed::spheres::planes::advance_mixed_sphere_scene`
automatically discovers sphere-pair and sphere/plane encounters in one native
mixed free/prescribed trajectory. `AcceleratingPlane` supplies a fixed normal
and a constant world material acceleration; its tangential velocity represents
conveyor motion. Every body pair and sphere/plane pair has an explicit material
law. The complete normal plane gap is quadratic, so the original ballistic
owner and interior extrema checks detect entries hidden between separated
endpoints. Body pairs retain the original accelerated swept-sphere owner.

All nearby features at an event enter one original mixed Coulomb solve. Plane
work uses the actual event-time material velocity, and native free recoil,
actuator transmission, prescribed support reaction and signed work remain in
the original physical accounts. Accepted endpoints also pass the original
body-pair path owner at zero duration, preventing a different floating norm or
rounded integrated pose from admitting a state that the next window rejects.
Bounded earlier-time retries change the candidate event duration; no gap,
position, velocity, normal or impulse is projected to pass admission.

`MixedSphereSceneResult::planes_after` rebases the original plane commands to
the accepted endpoint for repeated native continuation. Stable typed contact
keys distinguish body pairs from plane contacts. The old pair-only entry uses
the same event loop and retains its physical semantics. Existing bounds on
systems, links, contact count, events, native steps and geometry queries remain;
additional endpoint-owner checks consume the query budget. Up to 21 planes are
admitted, without silently dropping contacts when an active-set cap is reached.

The underlying enclosing carrier-sphere, coincident link-origin/COM, fixed or
revolute joint, and admitted centre-command restrictions still apply. Results
are floating analytic observations, not interval-certified CCD. Rotating
planes, arbitrary offset articulated surfaces, persistent contact, Zeno repair
and thermal loss partition remain outside this API. The focused tests in
`tests/articulated_mixed_sphere_planes.rs` cover analytic acceleration/conveyor
work, hidden entry and grazing, mixed simultaneous contact, covariance and
permutation, repeated continuation, resource/refusal atomicity and pair-only
compatibility.


## Sustained articulated plane forces through native motion

`articulated_contacts::step::sustained::step_sustained_articulated_planes`
couples the original force/displacement contact owner to original free and
prescribed spatial-motion owners. Every native trial starts from the same
accepted state and applies gravity, actuator effort, body wrenches and contact
point forces together at both midpoint stages. A resting force equilibrium
therefore balances gravity without a velocity reset or force-kick impact loss.
All points on each articulated tree retain the original full ABA coupling.

Planes have fixed frames and constant material velocity, including conveyor
motion. One held-force solve couples normal configuration complementarity and
Coulomb tangential displacement. Original relative material displacements,
plane source work, frictional dissipation, gravity and actuator/body work,
prescribed support work, actual momentum increments and independent source
impulses remain separate. The result distinguishes native integration error
from the defect between midpoint contact work and held-force displacement
work. No numerical defect is assigned to thermal heat.

The result is explicitly `ArticulatedSustainedStepCandidate<P>`: a physically
audited local candidate requiring caller-owned complete-path, omitted-feature,
contact-mode and any additional determinacy admission. Native midpoint/end
observations alone do not prove these properties. Strict final nonpenetration,
nonnegative normal force, cone/complementarity/passivity and work/momentum
checks apply; approaching impact states require an impact owner. A tiny signed
inactive remainder may trigger one exact-zero inactive warm-start restart of
the original solver, followed by fresh native replay and unchanged admission.
The original response/outer-iteration budgets are shared across both attempts;
no solved force, pose or velocity is clipped into acceptance.

Input states remain borrowed on all errors. Bounds retain 128 systems, 4096
total links and 21 contacts, with original native joint/state/source limits and
cancellation. This entry adds neither generic contact discovery nor a complete
Moreau--Jean scheme, a global Painleve/uniqueness classifier, frictional thermal
partition, rotating contact frames or interacting sliding body-pair surfaces.
`tests/articulated_sustained_contact.rs` covers exact resting support, moving
belts, raising planes, accelerating prescribed supports, zero-force detachment
and atomic path/impact/budget/cancellation/native-limit refusals.

## Continuous closed mechanisms on prescribed moving supports

`articulated::motion::spatial_loads::implicit::distance::prescribed::motion`
continues distance and spherical material closures together with affine joint
relations on the original prescribed base path. The initial constraint impulse,
midpoint joint rates and endpoint position closure solve simultaneously. A
second native impulse imposes endpoint tangency without changing position.
The original ABA, FK, normalized row mobility and nonlinear driver are reused;
only the at-most-32-row mobility is dense, with at most 4096 native links.
A representable endpoint-position chart defines the secant midpoint rates for
all trials and actual-impulse replay. Locked coordinates stay exact without a
final velocity projection; physical endpoint updates retain strict native
representability checks.

Diagonal and full passive six-coordinate bushings retain their original law.
The forcing entries resolve environmental link wrenches and actuator efforts at
the actual implicit stage, including the prescribed pose, twist and acceleration.
Callbacks receive elapsed h/2 and the absolute-clock state, with output buffers
reset from held loads. They must be smooth and repeatable on off-manifold probes
and stage their own history until success. Independent verification replays the
actual summed initial impulse and retains its accepted force sample through
native finalization. The open supported midpoint also exposes this same forcing
entry; no second force or support model is introduced.

`PrescribedBilateralLoopStep` retains both complete native impulse responses,
the smooth RNEA support receipt and the accepted midpoint. Total support impulse
and signed support work are compensated sums of these three independently
computed contributions. Constraint-event work, actuator/environmental work,
passive storage and damping remain separate. Actual per-link momentum and
kinetic increments are summed without subtracting large endpoint totals.
The smooth prescribed-path momentum and mechanical work defects remain visible;
no support reaction is inferred from a balance residual or converted into heat.

The caller owns a mutually consistent initial/midpoint/endpoint support path and
must resolve its time variation and nonlinear geometry. Initial/endpoint
position, endpoint tangency, native physical limits, row rank, nonlinear and
evaluation budgets, representable clocks and cancellation are admitted before
return. Borrowed inputs remain unchanged on every refusal. This adds continuous
supported closed motion, not contact activation, general branch certification,
or exact finite-step nonlinear energy/momentum conservation.

## Full passive bushings in rigid bilateral-loop motion

`implicit::distance::step_free_floating_bilateral_loops_with_coupled_bushings`
and its `_and_effort` entry point compose the existing full passive six-coordinate
stiffness/damping law with distance and material-point coincidence closures.
The original loop engine now consumes the existing `BushingLaw`; there is no
second constitutive or dynamics implementation. Cross-axis and
translation/rotation terms enter every simultaneous midpoint residual and the
independent actual-impulse verification. Existing diagonal and empty-bushing
paths retain their semantics.

Original reciprocal material-frame wrenches, full link inertias, carrier recoil,
both native constraint impulses, and separate actuator/storage/damping/work
accounts remain authoritative. Feedback uses the same repeatable trial contract.
The original matrix, geometry/chart, physical-limit, mobility, nonlinear-budget
and cancellation refusals apply. Finite-step constraint work remains visible;
this integration adds no energy-exactness, arbitrary active-law passivity or
general nonlinear-branch claim.

Four native test groups pass using an independent two-coordinate midpoint solution with
translation/rotation stiffness and damping, both carrier recoils and rigid-loop
reactions. Additional cases check diagonal/empty parity and atomic refusal.
All 21 existing distance, point-loop, bilateral-effort and coupled-effort cases
also pass with the original native owners. Full Cargo/DSR qualification and the
broader multibody bead remain outstanding.

## Joint gears and material loops in one free-floating solve

`articulated::motion::spatial_loads::implicit::distance::joint` composes the
existing constant affine `JointConstraint` rows with distance and spherical
material-point closures. Its diagonal/full-matrix bushing entry points, each
with optional repeatable joint feedback, use the same original midpoint and
native impulse owners. A single row mobility retains body/body, body/joint and
joint/joint responses, including complete link inertias and carrier recoil.
The initial impulses and midpoint rates solve together; one final native
tangent event enforces all rows together.

The combined bound is 32 scalar rows: a distance consumes one, a spherical
closure three, and a joint relation one. Only this bounded row matrix is
factored. Joint covectors enter the original native inverse-mass action; no
dense generalized mass matrix, separate trajectory solve or artificial mass is
introduced. Body rows retain their free-endpoint regularity test, and the
combined normalized reciprocal mobility must have resolved positive pivots.
A gear redundant with a geometric closure refuses instead of acquiring an
arbitrary reaction.

Initial and accepted endpoint joint positions/rates satisfy each row's original
Angle/Length units and absolute tolerances. Row coefficients scale against the
native angular/prismatic rate scales, so a consistent row rescaling also
rescales its nonlinear coordinates. Numerical trial states may leave the
manifold. The actual summed body and joint impulse events retain native
representability and physical state admission; passive reactions do not consume
motor effort ratings, while actual caller commands retain those limits.
Cancellation preserves the original loop error semantics.

`MixedBilateralLoopStep.bilateral` retains body-loop vectors in their existing
metre/N s coordinates, and both complete native impulse responses. Additional
fields expose joint multipliers, compact joint impulses, before/after row
residuals, and signed joint event work paired with actual mean pre/post rates.
That work is already included in the original total constraint impulse work,
including body/joint cross terms. A finite admitted tangent defect can do real
signed work; it is neither discarded nor counted again as damping heat.
World linear and orbital-plus-spin angular momentum remain independently
audited against external smooth impulses and actual reaction resultants.

Six focused native groups exercise an independent coupled row mobility with a
nonzero body/joint cross block, finite carrier recoil, both event-work
accounts and the joint-only boundary; a metre/radian geared mode with full passive matrices and feedback;
finite tangency work; a spherical closure plus angular lock against an
independently assembled full inertia tensor; empty-row parity, cross-family
rank refusal and cancellation; and passive reactions beyond motor ratings
with actual over-limit commands refused. All 25 existing affected loop/effort
cases also pass.

This is the original RATTLE-style bounded composition with measured finite-step
work, distinct from the continuous affine-reaction midpoint lane. It adds no
constraint activation, changing ratios, backlash, nonholonomic rolling,
arbitrary redundancy resolution, continuous-path enclosure, symplecticity or
energy-exactness claim. Full Cargo/DSR qualification and the broader multibody
bead remain outstanding.

### Retained spatial-gravity endpoint loads for coupled consumers

`gravity::step_in_field_with_loads` returns the original native step together
with the two actual conservative `GravityLoad` samples and the native state
after free drift, before the endpoint half-kick. `step_in_field` delegates to
that same implementation and keeps its original return type and map. Both
loads exclude the caller's separately held wrench; the native step includes it.

Coupled owners can therefore reconstruct actual half-kick work and external
linear/angular impulses and compare them with potential decrease without a
third source query or an inferred intermediate momentum. The original full
inertia, field/frame/domain validation and cancellation checks remain binding.
The receipt adds no field law, accepted-path domain guarantee or energy-exact
integration claim.

## Native impulse events on closed geared mechanisms

`implicit::distance::joint::apply_impulses` applies declared joint and body
impulses while existing distance, spherical and affine gear/lock rows react
together. It reuses the original objective geometry, combined body/joint
mobility, normalized Cholesky and native impulse application. The source
inverse-mass probe receives no physical endpoint speed admission; one summed
source/reaction event receives the original state and speed checks.

Positions remain unchanged and initial/final row residuals must satisfy their
original physical units and tolerances. Zero source still admits geometry,
tangency and combined rank, then retains the input motion with zero reactions;
an admitted finite defect does not activate a projection. Empty constraints
delegate to the original native impulse response. The bound remains 32 scalar
rows and 4096 links, without a dense generalized mass matrix or a contact law.

`BilateralImpulseResponse` separates constraint-only body and joint impulses,
material and joint multipliers, actual before/after residuals, and complete
native momentum/energy response. Actual source and constraint work use
compensated mean-velocity products. Objective row work and the total-work
partition residual remain separate observations, retaining rounding in the
wrench lift and summed impulses. No finite correction work is inferred heat.
The normalized reciprocal infinity-norm condition is estimated with retained
Cholesky unit-RHS solves for event consumers; continuous stepping does not pay
for those extra solves. Every inverse column and row norm is admitted finite.

Five focused native regression groups are added for independent geared
mobility/work and frame covariance; actual constrained endpoint speed versus
an inadmissible free probe; finite tangent work and passive exactness; spherical
and angular reactions against an independently assembled full tensor at an
admitted finite gap; and empty/rank/cancellation/retry behavior. Source/API and
analytical-oracle review completed. Native compilation, execution and rustfmt
remain unrun for this increment because the execution container was replaced
and the replacement has no Rust, DSR or RCH toolchain. These tests are not a
passing native or full-workspace qualification claim.

## Sustained reciprocal contact between finite articulated bodies

`articulated_contacts::step::sustained::step_sustained_articulated_contacts`
extends the same bounded original force/displacement solve and native
free/prescribed midpoint responses to finite body pairs and translating planes
in one contact inventory. The plane-only entry delegates to this shared engine.
Each finite pair applies equal and opposite world forces at its original link
material points. Native inertia, base recoil, joint coupling and prescribed
support reactions remain inside the original mechanical response.

Per-side spatial impulses and actual material-displacement work are retained.
Pair impulses are internal: actual gravity/body sources are rebuilt through the
original native spatial reducer at the native midpoint, before rounding, and
prescribed support and plane exchanges remain external. Relative witness
positions audit each pair's angular reciprocity before large world moments
can cancel. Numerical integration and contact quadrature defects remain
separate from physical frictional loss.

Initial, native midpoint and endpoint observations admit finite-pair normal-line
geometry and bound the actual force couple. Complete-path separation, frame
applicability, contact mode and additional determinacy evidence remain the
caller's responsibility. Fixed-witness sticking and central normal contact are
supported. Sliding surfaces needing newly selected material points, rotating
contact frames, generic collision discovery and a global Moreau--Jean theorem
are outside this local candidate's scope. No compensating couple, fictitious
inertia, state/force projection or automatic thermal credit is introduced.

G1/G3 regressions in `tests/articulated_sustained_pairs.rs` cover reciprocal
finite-body pushing, system/contact permutations, a moving three-body column
with either a prescribed support or a plane, guided finite-body sticking,
detachment, and artificial-couple refusals. The original sustained-plane
regressions remain tests of the same shared engine and its native refusal paths.


### Payload capture on closed geared mechanisms

`inventory::bilateral::capture_payload` joins the original rigid inventory owner
to the native mixed-constraint impulse event. The incoming payload's momentum
mismatch is formed against the actual carried velocity in one common
translating frame; the loaded full spatial inertia and all distance, spherical,
and ideal joint rows respond together. Only the combined physical endpoint
receives speed admission. Exactly carried payloads supply zero mismatch and do
not activate projection of an admitted finite tangent defect. Empty material
rows reproduce the original capture path exactly.

The inventory keeps its attachment identity, independent actual world and
projected momentum accounts, and nonnegative relative-motion capture quadratic.
Constraint-only physical wrenches are projected by the native zero-rate RNEA,
including any rounded base resultant. Actual common-frame constraint work is
reported and capped separately; it is not subtracted from capture loss or
inferred as heat. The mixed normalized reciprocal condition also observes the
existing payload condition floor, independently of pivot admission.

Five regression groups cover independently reduced geared-slider capture and
subsequent loaded motion, physical final-speed admission, world/frame and large
boost covariance, spherical closure with a full-tensor rigid-body reference,
finite constraint work and zero-mismatch passivity, exact legacy delegation,
and atomic condition/rank/energy/speed/cancellation refusals. Source/API and
independent mathematical review were completed. Native compilation, execution
of these regressions, and pinned formatting are **UNRUN** in this increment:
executor replacement removed the Rust toolchain and cached native environment.


## Automatic contact on translating articulated links

`articulated_contacts::trajectory::mixed::spheres::translating::advance_translating_sphere_scene`
accepts explicit `SceneSphere::TranslatingLink` surfaces alongside the existing
balanced `Carrier` geometry. Link surfaces attach to original native
material-point frames; their radii do not infer enclosure of link inertia or
omitted shapes. Each system uses one geometry mode. At most one sphere belongs
to a connected component of fixed joints, and every declared sphere pair and
sphere/plane pair has an explicit law, including distinct moving links within
the same tree.

Translating links require a prescribed base with exactly zero angular velocity
and angular acceleration, constant world translational acceleration, and
fixed/prismatic joints. Arbitrary fixed frame offsets/rotations, COM offsets,
full inertias, serial or branched topology, held native joint efforts/body
wrenches and world gravity retain their original meanings. Fixed link
orientations make the prismatic mass matrix constant. Native ABA and the
original point-motion map supply complete quadratic centre paths with genuine
cross-link mass coupling and finite recoil.

Before discovery, DD accumulation of original fixed joint axes times native
generalized accelerations, plus base acceleration, independently checks the
Cartesian acceleration. `UnresolvedLinkPath` refuses unresolved gravity/support
cancellation at the physical acceleration scale; it never substitutes a
repaired acceleration. Native joint admission also bounds interior coordinate
extrema: for a turning time t* in (0,h), q(t*)=q0+v0*t*/2 lies between q0 and
the existing midpoint predictor q0+v0*h/2. Initial/predictor/endpoint admission
therefore suffices without additional synthetic FK passes.

The shared original sphere/ballistic owners check complete paths. Actual native
endpoints and regenerated material witnesses retain strict nonpenetration and
bounded earlier-time retries. Simultaneous native contact retains joint/base
coupling, six-component support impulses, signed actuator/plane/base work, and
actual DD momentum/energy accounts. Limits, cancellation or late failure
publish no partial state. Existing bounds remain: seven spheres, 21 body
pairs, 21 planes, native link limits, 64 events, and explicit smooth-attempt
and geometry-query budgets; active-contact overflow refuses.

Coverage is floating analytic and applies only to configured spherical
surfaces. Rotating offset/revolute/helical link paths, omitted shapes, arbitrary
support commands, sustained contact, interval-certified CCD and thermal
allocation remain outside this entry.

Seven groups in `tests/articulated_translating_sphere_scene.rs` passed using
the original numerical owners in a bounded native Rust test graph: offset
drive/conveyor reactions, nonorthogonal serial mass coupling, same-tree
contact, mixed free recoil/covariance, hidden quadratic entry, atomic
joint-limit/cancellation refusal, and a finite-data gravity-cancellation
no-miss adversary. The twelve existing mixed sphere and sphere/plane groups
also passed. Formatting of all four changed source/test files passed. This
validation does not claim a full-workspace or DSR run.


### Retained-inertia release and ejection through closed mechanisms

`inventory::bilateral::release::release_payload` completes the same inventory
lifecycle for a free mechanism carrying distance, spherical and ideal joint
rows. The original attachment restores the exact retained link inertia; one
native mixed impulse event transmits the retained side's separation impulse
and all simultaneous row reactions. The separate payload responds through its
original full spatial inertia, with the same physical point, offset moments,
and free couple. Only the combined physical endpoint receives speed admission.

Passive and explicit zero-impulse release still admit geometry, tangent
residuals, combined rank and the existing condition floor, then preserve the
retained state bits. Active release keeps the actual common-frame separation
pair work and constraint work separately visible and independently capped.
The balance remains delta_K - W_separation - W_constraint, with no invented
capture loss, heat or fracture law. Constraint covectors use the shared native
RNEA projection. Empty material rows delegate to the original release exactly.
The auxiliary native event reports only retained-tree lab-frame source work;
the original release receipt accounts for both sides in its common frame.

Four regression groups cover an actual capture/loaded-motion/ejection cycle,
independent retained-mode recoil and final-speed admission, passive finite-slip
and signed-zero preservation, exact legacy delegation, spherical closure with
an independently assembled full-tensor and cofactor-inverse oracle, and signed
finite constraint work with separate budgets and atomic cancellation. Source/API
and mathematical review were completed, including an exact-rational check of
the geared-slider reaction and work equations. Native compilation, test
execution, and pinned formatting are **UNRUN** after executor replacement
removed the Rust toolchain and cached native environment.

## Simultaneous finite payload capture on a moving articulated tree

`articulated::motion::inventory::capture_payloads` accepts a nonempty set of
independent full-inertia payloads attached to any existing links at their actual
poses. Every attachment is measured against the original accepted configuration.
One cloned link inventory accumulates all loaded tensors, followed by one native
model construction and one free-floating articulated-body momentum solve.
Topology, link frames, joints, base pose and joint coordinates remain unchanged.
There is no sequential capture dynamics or dense generalized mass assembly.

Original per-link momentum and every incoming spin/orbital contribution use the
existing compensated spatial sum in one translating inertial frame. The final
loaded model supplies the coupled base and joint recoil. Existing ideal constant
gear/lock rows reuse their original prepared response on that loaded inertia
before final speed admission. A transient velocity from an isolated or partially
loaded capture cannot reject an admissible simultaneous endpoint.

The same native accounting reconstructs complete world and projected base/joint
momenta independently. Removed relative energy sums original retained-body and
incoming-payload velocity-change quadratics against the ONE final state, retaining
shared recoil terms. Common-frame endpoint energy and finite gear work remain
separate checks with independent allowances; neither rounding defect nor gear
work is assigned to heat. Singleton `capture_payload` uses this same owner.

Repeated target links return existing opaque release handles in input order.
Their exact prefix inertias describe inventory identities, not intermediate
physical motion. Handles on a given link release in reverse input order through
the original release owner; different links can release independently. The
physical result is permutation invariant to numerical accuracy, while handle
ordering follows the authored inventory order. Each combined mass increment
must remain representable under the original inertia admission.

Outside existing optional constraint response, assembly, momentum and energy
passes require O(links + payloads) work/storage. Additional loops poll every 64
entries, and indivisible native tree passes retain the 4096-link bound. All
failure paths leave inputs unchanged and publish no partial model, motion or
attachment. Contact discovery, overlap validity, finite grasp dynamics, capture
of an independently articulated mechanism and thermal conversion remain outside
this instantaneous rigid attachment model.

Five focused G0/G1/G3/G4 regression groups cover analytic different-link recoil
and small loss under a billion-scale common velocity; opposite incoming spins
whose individual captures violate a speed bound but joint capture satisfies it;
reverse repeated-link release recovering exact retained inventory; independent
COM/screw and all-piece parallel-axis sums for mixed joints and full tensors;
input permutations and fixed-scene/boost covariance; singleton equality and
existing locked-joint response; and late inventory, work and final-cancellation
refusal with exact retry. Pinned nightly Rust formatting/parsing completed.
Native compilation and numerical execution remain unrun in this session.


### Simultaneous payload capture through closed mechanisms

`inventory::bilateral::capture_payloads` connects the original many-payload
inventory owner to the mixed distance/spherical/gear impulse event. Every
incoming-minus-actual-carried momentum mismatch is accumulated by the existing
compensated spatial sum, then all loaded inertias and constraint cross responses
enter one native event. Singleton and batch capture share this lifecycle; the
obsolete separate singleton lifecycle is removed. Empty material loops retain
the original batch path.

The same final state determines every outgoing attached velocity and the full
nonnegative capture-loss quadratic. Opposing incoming mismatches can cancel
the mechanism's source while still removing relative energy; such an event
preserves admitted finite tangency defects. Exact release-handle stacks,
common-frame actual constraint work, condition policy and endpoint speed
admission retain their original owners.

Four focused regression groups cover independent two-payload reduction and
loaded continuation, combined endpoint limits, repeated-link inventory with
opposing arrivals, large common boosts and finite constraint work, exact
single/empty delegation and atomic refusal/cancellation. Source/API and
independent exact-rational oracle checks were completed. Edition-2024 parsing
and formatting pass with available nightly-2026-08-31 rustfmt; the repository's
pinned nightly-2026-07-06 is unavailable. Native compilation and numerical test
execution remain **UNRUN**: no existing FrankenSim native dependency graph,
DSR or RCH is available in this runtime. This is not a native test-pass claim.

## Finite rigid-component transfer networks

`assembly::transfer::rigid::transfer_rigid_components` redistributes original
finite components between moving rigid assemblies in one atomic event. Each
leaf has exactly one destination; every final body remains nonempty. Cyclic
handoffs use only original inventories, including complete inventory exchange.
The final body axes remain the corresponding original receiver axes.

A donor's components assigned to one receiver remain a rigid group. The native
`AssemblyPartition` prepares full group inertias, and its original release
owner applies all authored force/couple pairs together. Shared-group recoil
cross terms therefore enter the source work once. Every pair must cross released
groups. After all donors release, the original `RigidAssembly::latch` admits
each receiver's complete incoming set, including intrinsic spin, orbital momentum,
and nonnegative relative-motion loss. No leaf-release/relatch surrogate is used.

Final inventories retain flat original component geometry and properties,
alongside original body/part indices, so later transfer events retain material
identity. Flattening changes bookkeeping only: native aggregate inertia and
momentum are kept. World component poses do not change during this instantaneous
event. Contact discovery, attachment strength and geometric admissibility,
finite flight, actuator energy debit and thermal deposition remain caller-owned.

One reusable destination map and donor/part-ordered groups give linear routing
and staging in bodies, components and impulse pairs; there is no dense network
mobility matrix. Aggregate bounds are checked before network staging. Single
groups skip release, singleton unchanged-axis receivers skip latching, and
isolated bodies preserve their exact original state and inventory.

Native per-phase checks and the existing polarized common-frame transfer audit
check actual final stored mass, linear/angular momentum and kinetic change minus
source work plus sticking loss. Source and loss accounts remain separate per
body and globally. Late refusal or final cancellation publishes no endpoint.
All five rigid-component-network integration tests now pass in the selected
native source build alongside the prepared local-transfer qualification below.
Full-workspace Cargo and distributed gates remain unrun.

Five focused G1/G3/G4 regressions cover independent shared-recoil work and
intrinsic-spin sticking loss, retained packets without artificial leaf loss,
full-tensor cycles against a flat-leaf mechanical oracle, scene/boost and body-
permutation covariance, exact no-op and whole-inventory reuse, and late invalid
inventory, aggregate caps and final cancellation/retry. The analytic fixture's
source work is 41/5 J and sticking loss is 117/56 J; rational arithmetic checks
those values independently of the native Rust implementation.


### Simultaneous release with shared retained-mechanism recoil

`inventory::release::release_payloads` and
`inventory::bilateral::release::release_payloads` detach a nonempty payload
collection in one physical event. Per-target handles validate exact retained
inventory prefixes in reverse capture order; this ordering creates no
intermediate dynamics. All outgoing bodies inherit the original loaded FK,
all selected inertias are restored before one native tree response, and each
payload receives its own full spatial-inertia response. The existing singleton
release entries now unwrap the same batch owner. Supported inventory continues
to reuse the original attachment and separation-wrench admission helpers.

Compensated per-target source sums retain shared recoil. Passive and exactly
cancelling tree sources preserve the retained state bits after constraint
admission, while individual nonzero payload impulses still act. Every source
pair uses the same final retained velocity, retaining cross terms and signed
per-pair work; its DD total is independently capped by the existing absolute
source-work allowance. Mixed material/joint reaction work and actual energy
and momentum defects remain separately admitted. No heat or fracture law is
inferred. Assembly/accounting outside the original bounded constraint solve is
O(links + payloads), without sequential model rebuilds or a dense mass matrix.

Five regression groups cover analytic shared recoil and a sequential-release
counterexample, input permutations, exact empty delegation, reverse inventory
and atomic cancellation, combined mixed-constraint endpoint limits, finite
slip and signed work budgets, cancelling sources with nonzero payload recoil,
and full-tensor/off-axis independent inverses with frame/boost covariance.
Source/API review and exact-rational scalar checks completed. The four changed
release source/test files pass edition-2024 parsing and formatting with available
nightly-2026-08-31 rustfmt. Pinned July formatting, native compilation and
numerical execution remain **UNRUN** because this runtime has no original
FrankenSim native build graph, DSR or RCH. No passing native qualification is
claimed.


### Native inverse-mass actions for articulated multiphysics composition

`articulated::free_floating_impulse_delta` and
`prescribed_base_impulse_delta` expose the existing O(links) zero-bias ABA
actions, with original `ImpulseVelocityChange` and `PrescribedImpulseDelta`
values. Their numerical implementations are unchanged. Body impulses are
about link origins in link axes; joint impulses retain their revolute/prismatic
units. The prescribed action returns its original condensed support impulse.

These are linear mobility probes, independent of incoming velocities and
smooth loads. They validate configuration, dimensions, finite inputs and
native inertia admission, while permitting probe velocities beyond physical
joint speed limits. Coupled consumers must submit their combined impulse to
the original free/prescribed physical impulse owner for actual endpoint and
momentum/work admission. The articulated/modal consumer exercises that
distinction with an admissible combined response and an atomic actual-speed
refusal, alongside independent coupled-mass and prescribed-support oracles.
No state is changed by a mobility probe and no finite force duration is inferred.

## Prepared local finite rigid payload transfer between moving bodies

`assembly::transfer::rigid::PreparedRigidPayloadTransfer` prepares a nonempty
proper subset of an original rigid assembly for local handoff to another finite
body. The existing rigid-component network API retains its general routing and
flat inventory interfaces. The prepared local entry accepts a native finite
receiver directly, caches the two outgoing groups, and admits explicit matching
material attachments. The existing partition owner retains complete original
geometry and full central tensors for the remaining donor and payload. Every
material component belongs to exactly one group, and the donor retains positive
mass throughout. Preparation is O(parts); an event is O(separation descriptors),
with two finite bodies in each native solve. There is no per-leaf reassembly, artificial
release/relatch loss within a retained group, or dense multibody system.

One original native partition release applies all equal/opposite force impulses
and free couples together. Original component indices and donor-COM application
points select the two retained groups. Their complete recoil includes shared
impulse cross terms in the separation work. One original native rigid latch then
captures the whole payload at its existing pose, preserving intrinsic spin,
orbital momentum and full inertia. The receiving assembly keeps the original
receiver axes and contains the original receiver and whole payload as its two
parts. The prepared partition retains payload leaf geometry and original indices
for later expansion or subsequent inventory events. No free payload state is
returned alongside the captured receiver.

Both attachment arms refer to the original body COMs and axes. Compensated
origin-relative coincidence admission precedes release; the actual released
payload's local and rounded world points must also pass before latching. Equal
rounded coordinates at a remote world origin cannot hide different local
attachments or an unrepresentable released COM. This is an instantaneous rigid
handoff, with contact discovery, overlap, connectivity, strength, finite grasp
dynamics, free flight and material identity remaining caller-owned.

The shared transfer envelope admits all original donor components plus one
receiver before event staging. Every separation descriptor, including a zero
one, consumes an impulse slot; the receiving latch consumes one further slot.
Native release and latch gates precede a combined check of actual final mass,
linear momentum, spin plus orbital momentum and common-frame kinetic balance.
The shared balance owner accounts for each actual output orientation, including
native quaternion normalization. Separation work and nonnegative sticking loss
remain separate: the residual is delta_K - W_separation + E_sticking. A caller
must debit the work source and select any thermal destination explicitly.
Cancellation or a late refusal publishes neither endpoint and changes no input.

Six regression groups independently verify finite spin recoil and closed-form
separation work/sticking loss; off-center multi-part payloads with independently
assembled full tensors and a cofactor inverse; world spin/orbital conservation
and all-leaf inventory recovery; rotated, translated and commonly boosted
scenes; proper subsets, positive retained mass and exact work-cap boundaries;
remote-origin geometry refusal; and every observed preparation/event
cancellation point through final publication, late balance refusal and retry.
The original point-transfer test's quaternion helper now uses the public
components API, allowing its existing covariance fixtures to compile.

Pinned nightly native compilation and all 37 selected integration fixtures
passed: six rigid-payload, seven assembly, six partition, seven point-transfer,
six point-transfer-network, and all five concurrent rigid-component-network
tests. The native harness uses the actual core and assembly source with
module-path adaptation and genuine cached dependencies;
no physical solver is replaced. Full-workspace Cargo and distributed gates
remain unrun in this session.

## Prescribed final-relative-twist payload ejection

`inventory::release::ejection::eject_payloads` commands each payload's final
angular and origin-linear velocity relative to the FINAL retained-link field
at that same payload origin, in the payload's original axes. This preserves the
difference between an origin velocity and a COM velocity. Original opaque
attachment chains and full spatial inertias remain authoritative; command
count is capped at 4096 before staging.

For payload map S_i and spatial inertia I_i, loaded momentum gives
M_loaded*delta_v = -sum(S_i^T*I_i*d_i) + constraint_impulse.
The existing native ABA and original row response determine shared recoil.
Derived payload impulse I_i*(S_i*delta_v+d_i) uses small velocity CHANGES,
never subtraction of large before/desired endpoint velocities. Force is placed
at the original payload origin with its spin impulse represented as a free
couple; this is a wrench decomposition, not inferred actuator geometry.
No dense 6N mobility or independent pairwise ejection sequence is introduced.

The loaded-side source uses the same native point/couple decomposition and
compensated link accumulation as physical release, retaining exact-cancelling
source semantics. True all-zero commands skip auxiliary recoil and preserve
passive native state bits, including admitted finite tangent defects. Active
rows retain their actual initial defect correction. The mixed distance/
spherical/gear entry `inventory::bilateral::release::eject_payloads` removes only
uniform translation from its auxiliary state, retains original angular/joint
rates, and recovers the native full-body delta from summed source/reaction
impulses. Auxiliary work is never counted as physical separation energy.

The original simultaneous retained-inertia release owns all returned states,
source/constraint work, momentum/energy and physical speed admission. Both
actual relative linear and angular endpoint residuals must then meet declared
component tolerances before publication. At an exact commanded endpoint, pair
work is 0.5*d_i.pairing(payload_impulse_i); cross terms from common recoil remain
inside that impulse. The original separate absolute work caps apply.
Outside original constraint handling, work/storage is O(links + payloads).

Focused regressions cover the independent 213/28 J shared-recoil oracle,
full-tensor off-center payload-origin commands and scene covariance, mixed
distance/gear recoil including finite admitted tangent defects, passive state
bits, and target/inventory/work/cancellation refusal with exact retry. Pinned
Rust parsing is separate from native compilation and numerical execution,
which remain unrun while local storage is full and DSR/RCH are unavailable.
Contact/fracture decisions, real actuator placement and finite-time forcing,
supported-base ejection, mechanical-source evolution and heat remain external.

### Supported mixed material/joint impulse events

The prescribed-base boundary now composes existing distance/spherical material
closures and constant gear/lock rows in one native event through
`motion::spatial_loads::implicit::distance::prescribed::apply_impulses`.
Native prescribed ABA supplies all material/joint cross-mobility columns and
one final summed source/reaction response. Fixed/driven supports and mobile
scalar root joints retain their actual geometry, full spatial inertia and
independent six-component support reaction. Only the final physical endpoint
receives native speed and representability admission.

The original material geometry, finite-gap spherical couple, normalized
reciprocity/rank admission and retained Cholesky factor are shared with the free
mechanism. No dense generalized mass matrix, surrogate support mass or impulse
to force conversion is introduced. Zero source admits existing geometry,
tangency and full mixed rank without correcting finite admitted rate defects;
empty rows delegate to the native prescribed impulse owner. The source and
constraint work partition excludes signed support work, which remains in the
native receipt. Rows do not engage themselves and no time is advanced.

Five focused regressions were authored in
`tests/articulated_prescribed_bilateral_impulses.rs`: independent scalar coupled
mass references, fixed/driven support and world covariance, mobile roots,
full COM/off-diagonal tensor and spherical finite-gap couples, finite-slip
reaction work, empty/zero-source behavior, joint limits, combined rank,
cancellation and deterministic retry. Source/API and independent analytical
review completed. Available nightly-2026-08-31 rustfmt parsed the changed parent
without rewriting it and formatted/checked both new files; this is not the
repository-pinned toolchain. Native compilation and numerical Rust execution
remain **UNRUN** here because nightly-2026-07-06, the complete dependency graph,
DSR and RCH are unavailable. No benchmark, all-runtime or full-plan completion
claim is made. The broader `frankensim-ext-mbd-core-vqqt` bead remains open.

### Supported payload capture and release with material closures

`motion::inventory::supported::bilateral` connects supported rigid payload
capture, passive release and point-impulse/free-couple separation to the native
prescribed mixed material/joint event. The existing inventory lifecycles remain
the owners: capture combines the original full rotated payload inertia before
the simultaneous mismatch/reaction response; release restores the exact retained
inertia before the simultaneous recoil/reaction response. The new path admits
distance/spherical closures and gear/lock rows together, including mobile root
joints, and retains the authoritative native support impulse and signed work.
The original unconstrained and affine-only public entry points retain their
existing branches.

Actual body and joint constraint impulses now share the original zero-rate RNEA
projection and common-frame work helper. The base covector uses the current root
joint transform; old fixed-root callers keep their original transform and f64
work boundary. Supported inventory subtracts the full projected joint impulse,
not only the affine A^T multiplier. The new reaction receipt exposes the full
base/joint covector, common accounting boost, actual midpoint work and final-rate
pairing, separately from the native auxiliary event's lab-frame work.

Capture retains its nonnegative physical relative-motion quadratic Q and the
identity delta_K - W_support + Q = J_constraint dot v_final, including the
independently bounded finite tangency defect. Midpoint reaction work is capped
separately and is not subtracted from Q. Release retains unrounded Dd constraint
work until subtracting support, source and reaction work from actual kinetic
change. The supported separation policy still refuses negative source work.
Passive/zero-impulse release admits geometry, tangency and full mixed mobility
on retained inertia while preserving the original state and applying no
constraint correction. Existing representability checks still refuse nonzero solved velocity changes
erased in actual endpoint body motion, including possible roundoff at a fully
immobilized driven configuration. No hidden rate clipping is introduced.

Five focused groups in `tests/articulated_supported_bilateral_payloads.rs`
cover an independently derived capture/passive/ejection cycle, simultaneous
endpoint limits, finite-slip work and caps, exact passive state, full-tensor
payload response with spherical closure and a nonzero mobile-root coordinate,
independent virtual-work projection, rotated/translated frames, driven support
work, a large common boost, combined-rank/condition/budget/token refusals and
cancellation. Root and independent review found no blocking source/API defect;
exact rational scalar and finite-slip work checks passed. All changed consumer
files and the new tests pass the available nightly-2026-08-31 rustfmt check.
Native compilation and numerical Rust tests remain **UNRUN** under the previously
recorded toolchain/dependency limitations. No throughput, full supported
continuous-loop integration, supported batch-inventory or full-plan completion
claim is made; `frankensim-ext-mbd-core-vqqt` remains open.

## Finite physical nozzle transport on full-tensor motion

`rocket::rigid::nozzles::transport` adds
`TensorRocketState::advance_nozzles_transport`. A private observer consumes the
original accepted rigid stages; native motion, mass/inertia, burnout and work
limits retain one owner. Held nozzle rates, geometry and complete local relative
energy belong to `NozzleTransportForcing`; additional spatial loads stay separate.
There is no emission after burnout or a second application of propulsion.

The original analytic velocity moments integrate physical linear export and COM
kinetic transport. Native midpoint spin supplies the additional co-rotation and
tangential-jet energy correction. Local relative energy includes jet kinetic and
the supplied caloric datum once; its integral and both frame corrections remain
separate. No engine, pressure-flow, internal-flow or thermal closure is inferred.

Two-point orbital quadrature follows the same analytic COM trajectory. It is
anchored at each native stage and translated with its analytic linear impulse;
the requested fixed-world origin is applied once using the integrated impulse.
Spin export remains the original native angular outflow integral. Returned
momentum residuals retain actual endpoint and quadrature defects. The complete
rotating/depleting composition is second order and requires temporal refinement.
Work/storage per native substep are constant, independent of nozzle count;
cancellation or exhausted work publishes neither a partial body nor stream.

Five focused native tests cover analytic forced exhaust through burnout, native
endpoint parity, an independent continuous rotating-particle reference with
second-order momentum/angular-momentum/energy refinement, origin/boost changes,
signed caloric data, dry coast, work limits and late cancellation. The bounded
actual-source native run passed all 69 included tests; workspace Cargo gates were
not run in this constrained environment.

## Finite actuator power on original rotor motion

`rotordynamics::transient::RotorMidpoint::step_power_limited` advances the
original full mass, gyroscopic, damping and stiffness operators under one
finite commanded generalized-effort pattern. It selects the largest admitted
fraction in `[0,1]`, retaining the command's force/torque proportions and never
exceeding an individual commanded magnitude. Separate nonnegative motoring and
regenerative limits constrain the signed net interval-average mechanical power.

The same cached nonsymmetric midpoint factorization supplies the unforced
response and command influence. Their actual source response is the scalar
quadratic `P(alpha) = alpha * (b + alpha * a)`. Both roots are considered when
necessary: active bearing coefficients and velocity reversal need not make
source work monotone in drive fraction. The pure-quadratic start-from-rest
case takes separate square roots before division, preserving a representable
small effort fraction even when the corresponding power ratio underflows.
Initial speed is never a divisor, so finite commands remain finite at stall.

The selected effort enters the original `step` once. Its rounded physical
endpoint, momentum admission, source work and independent kinetic, bearing,
stiffness and gyroscopic work remain authoritative. Final actual source work
per second must meet the signed limits within the explicitly authored
`absolute_power_tolerance_w`; the original work/energy residual per second must
also meet that allowance. Refusal or cancellation returns no partial state or
source debit. All inputs remain borrowed and unchanged.

Preparation and original state ownership are unchanged. The additional work is
two bounded cached linear solves and fixed-size stack arrays, followed by the
original accepted solve; each step stays `O(n^2)` for at most 32 coordinates,
with no refactorization, matrix inverse or nonlinear search. This is one net
mechanical source, without per-coordinate power allocation, instantaneous
power limits, finite source storage, electrical efficiency or thermal losses.
An interval may absorb and return energy while its net work is zero; within-step
reversal and source constraints require a resolved time grid and their own laws.
The original constant-spin linear rotor model and temporal-error limits remain.

Six focused native tests cover finite stall torque, unchanged full-command
parity, constant-power acceleration and regeneration with exact energy and
three-rung position refinement, nonmonotone/active-bearing root selection,
full tensor/gyro/cross-coupling and coordinate-unit changes, `1e-180..1e180`
power scaling, a representable `1e-300` effort fraction, strict final energy
refusal, and every observed cancellation boundary. All six pass alongside all
12 original rotor operating-point tests and four unchanged transient
regressions in a source-selected native build using the original numerical
owners. Unrelated rigid-disk/periodic module exports were omitted from that
build; whole-crate Cargo/DSR verification is not claimed.

### Peak actuator power during acceleration and reversal

`RotorMidpoint::step_peak_power_limited` separately constrains both signed
endpoint powers through `RotorPeakPowerControl`. Under held effort and the
midpoint method's linear velocity reconstruction, power is affine in time;
these two endpoint limits therefore admit every reconstructed intermediate
power, including a reversal whose net work is zero. The ordinary average-power
API retains its original semantics.

Initial power first bounds the effort fraction to one interval. Final power
is an original-response quadratic on that interval; the stable existing
two-root solver chooses the largest feasible fraction, even when active
bearings create disconnected feasible sets. Native full M/G/C/K, finite stall
torque, cached factors and O(n²) stepping are retained. Actual endpoint powers
are evaluated independently after the accepted native step, with signed
source work and the original work residual; refusal publishes no motion.

The peak guarantee concerns the declared time interpolation, not unresolved
extrema of the exact continuous trajectory. Time refinement, source storage,
electrical efficiency and thermal losses remain separate physical owners.
Focused regressions cover stall acceleration, zero-net-work reversal,
regenerative limits, active bearings, original tensor/gyroscopic coupling,
coordinate congruence, extreme power scales and atomic cancellation.

## Moving articulated trajectories with coupled dry friction

`rigid_rattle::distance::graph::joint_friction::relaxation::trajectory::advance_joint_friction_trajectory` advances actual native poses and momenta by composing the original constrained graph half steps with the continuous coupled-friction phase. Geometry and the complete native constrained mobility are rebuilt at every moved midpoint. The original flow travel remains a path-work measure; it is never applied as a second pose drift.

A fixed stack source bank permits changing authored capacities and drives at the actual midpoint states and physical time, preserving joint indices and kinds. No bearing normal pressure is inferred from aggregate bilateral reaction. Strictly ordered caller-authored breakpoints split known load switches and short pulses; unknown features still require a suitable maximum step. The represented endpoint-minus-start duration owns both clock advancement and integrated force/work. Interval partitioning integrates a small terminal remainder instead of discarding it.

One-split versus two-half-split differences control positions, orientation, world spin and linear momenta, per-port scaled travel, actions and work. Raw differences use no assumed-order Richardson divisor; stopping can reduce splitting order. The accepted two-half paths accumulate numerical estimates against caller budgets, not certified global error bounds. Trial/step/time/workspace limits and representability failures refuse without partial publication.

The final receipt retains original mechanical energy including uniform-gravity potential, original free-step force/displacement and body-torque/midpoint work, physical axial drive work, nonnegative friction loss, signed constraint work and its absolute accepted-step sum, and the raw balance residual. These independent accounts remain separately gated. Exact static balance belongs in the axial drive ports: splitting separately authored body wrenches can require further refinement or refuse at the configured limit. There is no inferred heat/wear partition or distinct static/kinetic coefficient in this owner.

Nine direct regressions cover finite-carrier stopping distance, reversal poses/work, analytical coupled-rotor attitudes, authored short pulses, pose-dependent spring drive, a rotating rail with rotated full body inertias against an independent centrifugal/friction ODE, gravity and external work, large-epoch represented time, and atomic late-source/cancellation/refinement refusals.

## Native joint breakaway and sliding friction

The continuous joint flow now accepts `RevoluteStiction` and
`PrismaticStiction` sources with independently authored static and kinetic
capacities. Their complete constrained mobility, bilateral reactions,
full-tensor inertia response and native momentum application remain with the
original graph owners. Legacy sources and explicitly equal bounds retain the
original common-bound path. The six-port distinct-bound path uses the shared
contact mode resolver, with a checked total bound of 729 mode evaluations per
allowed affine segment.

Static and kinetic force/torque conversion uses the same work-conjugate scale as
drive. In particular, exact authored static balance is not broken by rounding
the static force downward alone. Conservative integrated-impulse admission and
actual native residual audits remain in force. Reports expose mode diagnostics
and complete path speed extrema including measured native endpoints. Coupled
mode ambiguity remains an explicit refusal. Underlying segment and impulse
budget refusals preserve their original error identity so adaptive trajectory
refinement continues to operate.

Existing motion source callbacks can now hold a joint, trigger breakaway, and
continue at the lower sliding force. The law recognizes exact rest and isolated
resolved stopping events; it does not silently capture a small nonzero native
velocity or supply persistent geometric locking or compliant presliding state.

Eight tests exercise static versus moving response, tiny real slip, exact
thresholds across nonbinary hinge scales, recapture and reversal, coupled rotor
reactions above the kinetic limit, source permutation, full-inertia eccentric
recoil and frame/boost covariance, exact legacy equivalence, existing adaptive
motion and segment-budget refinement, invalid sources, budgets and cancellation.

## Ordered dry-interface cards on native joint flow

The continuous flow's `interface` child binds original `DryInterfaceSystemCard`
objects to ordered native joint bodies and explicit physical load geometry.
Uniform axial sliding derives capacities as mu*N. A fully sliding cylindrical
sleeve derives torque as mu*N*r and surface speed as r*abs(omega), using an
authored physical radius independent of the hinge's numerical angular length.
Normal force, nominal area, temperature and radius are held source data;
pressure is explicitly nominal N/A, not an inferred distribution or bilateral
reaction magnitude.

The constant Coulomb rung supplies distinct static and kinetic coefficients.
Other constitutive laws refuse rather than being replaced by zero-speed values.
The exact borrowed ordered card, surface history, provenance and authority are
retained in the receipt; this does not promote caller data to fs-matdb admission
or empirical validation. Original temperature and pressure applicability is
checked, then complete path speed extrema and direct native endpoints are
validated before returning the candidate. Both lower and upper speed bounds
apply, including arrests and intermediate peaks. Late domain refusal or
cancellation publishes no partial body state.

Six tests cover rest versus slide, physical radius and conditioning invariance,
pressure/temperature/lower-speed refusals, an internal speed peak with safe
outer endpoints, binding permutation and ordered card identity, and unsupported
laws/geometry/cancellation. The shared kernel ambiguity test additionally checks
incomparable maximal sticking sets of different cardinalities.

### Uniform gravity work in free-motion receipts

`FreeFloatingMotionStep::gravity_work` exposes the uniform-gravity contribution
already computed by the explicit and implicit native motion owners. It is
included in `external_work`, and must not be added to it again. Both owners
retain their existing compensated mass-times-gravity-times-COM-displacement
calculation before adding the world origin. Conservative coupled storage gates
can use its negative as the body gravitational-potential change without
subtracting large absolute pose potentials. The integration equations, force
stages and momentum reconstruction are unchanged.

## Adhesive rolling spring on a native rigid pair

`tangential_friction::adhesive::rolling` connects an original
`SourceJkrRollingContact` to the actual finite-body inverse inertias along the
contact frame's longitudinal bending axis. An analytic elastic midpoint
predictor is accepted only after one unchanged source evaluation confirms its
residual and absence of yield. Other cases use a bounded monotone scalar solve.
Both paths couple the source's exact work-mean moment to a discrete-gradient midpoint
impulse phase. Both bodies receive equal/opposite free angular impulses through
the existing native resultant path. No normal or tangential force is invented.

Sealed history binds source, pair indices, material arms, frame, poses and mass
properties. A changed geometry cannot inherit old elastic rolling storage.
Source evaluation and native preparation never mutate accepted state; the
publication helper commits both momenta and rolling memory after cancellation,
angle, velocity, individual/pair momentum and work checks. Compensated actual
kinetic change retains weak rolling work beside large common translation.
Recoverable rolling energy and irreversible yield loss remain separate, and
neither includes the original normal JKR storage or an assumed heat partition.

The constitutive increment is a monotone scalar path. The inertial update is a
time discretization. Endpoint reversals require a zero-loss elastic trial and
a guarded initial mode-energy bound strictly below yield storage; this excludes
hidden yield anywhere on the elastic orbit. Absolute roundoff allowances use
original angular momenta and minimum principal inertias to raise the rate and
lower mobility, retaining cancellation and conditioning effects. The fixed
4096-epsilon margin is independent of authored audit tolerances. Unresolved,
subnormal or nonfinite bounds refuse; this is a conservative floating-point
admission screen, not an interval or time-accuracy certificate. Other reversals
still require refinement. This phase does not enforce no-slip rolling, integrate
poses, evolve normal area, transport memory or supply general vector rolling. G1/G3/G4
cases in `tests/adhesive_rolling.rs` cover analytic elastic pair recoil, tensile
yielding, full inertia in rotated frames, weak work beside common translation,
and atomic cancellation, geometry and numerical-budget refusals.

## Prescribed normal loading with retained adhesive rolling memory

`tangential_friction::adhesive::rolling::normal` extends the native rolling
phase with a prescribed attached JKR indentation update. The original normal
source determines the new contact area while the exact elastic fraction and
accumulated rolling angle are held. Material, curvature, critical rolling shift,
calibration and small-shift applicability remain bound to the original source.
The resulting radius-dependent stiffness then drives the existing native
midpoint recoil using both actual body inertias.

Receipts separate original normal-contact work, signed rolling parameter work,
applied body-wrench work, recoverable storage and irreversible rolling loss.
Changing stiffness at held elastic shift supplies or returns rolling energy; it
does not manufacture dissipation. Original factored incremental work/storage
primitives retain weak changes beside large preload. Independently observed
absolute storage differences are cross-checked with an explicit rounding
allowance, while native mechanical work retains its original acceptance gate.

Preparation returns an immutable candidate. Publication updates both native
body momenta, the original-source continuation and its sealed rolling history
together, only after all domain, geometry, work and cancellation checks.
Below-fold indentation, changed geometry, invalid wrenches or cancellation
cannot publish a partial area or history update.

This is a prescribed normal stage followed by rolling recoil at held body poses
and frame. It does not solve normal body motion or its impulse, rolling-induced
changes of equilibrium contact area, detachment, frame transport, material
refresh or caloric deposition. The tests in `tests/adhesive_rolling/normal.rs`
exercise analytic area-dependent recoil with retained preload, exact replay of
an unchanged normal state, and atomic refusal/retry of all participating owners.

## Bonded adhesive torsion on a native rigid pair

`tangential_friction::adhesive::torsion` binds the original JKR torsional source
to a finite rigid pair and the original contact normal. Both actual inertia
tensors determine normal-axis mobility; equal/opposite free couples change
both original angular momenta. Original twist storage and irreversible loss
close the compensated native kinetic change, including under large common
translation. The source's explicitly calibrated critical twist and mode-III
release/adhesion approximation ceiling remain load-bearing.

Rolling and torsion share one private midpoint phase: the source-verified
analytic elastic predictor, bounded monotone solve, guarded elastic-turn
screen, and independent actual momentum/work audits. Their public source,
history, axes and reports remain distinct. History binds pair, arms, poses,
frame and mass properties; cancellation and any failed solve/source/audit leave
bodies and memory unchanged. This fixed-normal phase does not enforce force
reactions, transport rotating contact history, evolve area, resolve local
partial slip, combine sliding capacities or assign irreversible loss to heat.
G1/G3/G4 tests cover elastic recoil with rotated full inertia, tensile yielding,
common-translation invariance, safe elastic turns, hidden-yield refusal and
atomic source/geometry/budget/cancellation refusal. Existing rolling tests
remain the regression seam for the shared native arithmetic.

### Applied loads on adhesive rolling and torsion

Both native moment owners provide `prepare_*_with_wrenches` and
`apply_*_with_wrenches`. The two original bodies' held world COM forces and
body-frame torques enter the same source-owned midpoint solve. Actual full
inertia determines the free angular acceleration; an applied scalar effective
torque or an already rounded free-body predictor does not replace it. Contact
and applied impulses are summed before updating each body once. Existing
zero-load entry points use this same path with zero wrenches.

Reports retain actual COM and held-frame angular travel, external work, contact
work, recoverable storage and irreversible rolling/twist loss separately.
Independent checks close kinetic change against both work channels and total
kinetic-plus-storage-plus-loss against external work. Individual commanded
impulses and pair linear/orbital-plus-spin changes are checked against the
actual applied loads. Applied force does not silently become a contact reaction.
Poses, normal area and contact history geometry remain held; returned angular
travel is not a finite-rotation pose update. A reversal under an applied body
torque requires interval refinement because the unforced energy-orbit bound
does not constrain a driven orbit.

The existing rolling/torsion G1/G3/G4 tests additionally cover analytic driven
elastic response with rotated full inertia, loading from rest into irreversible
yield, independent applied work and angular momentum, common translation,
zero-load equivalence and atomic invalid-load/cancellation/reversal refusal.


## Finite battery power in original rotor motion

`rotordynamics::transient::drive::voltage::battery` connects the original finite battery/converter to `RotorVoltageMotor`. Its cached mechanical Schur response gives an affine midpoint winding current at the converter's discrete output voltage. That current and its exact tangent enter the original six-state electrical solve; the accepted voltage recovers the original full-coordinate M/C/K/G motion. There is no second inertia, staggered voltage/current interface iteration, clipping after acceptance, or energy correction.

`BatteryRotorStep` retains the original source and mechanical records, checks current and work continuity, and closes electrical/magnetic/kinetic changes plus original C/K/gyroscopic work, electrical/winding heat and auxiliary export against external generalized-load work. The nine heat channels are the original eight battery/converter channels followed by winding Joule heat. Active damping and nonsymmetric stiffness retain signed operator work and are not relabeled heat. The direct winding connection requires unit empirical bridge efficiencies because converter losses already have physical owners.

Voltage curtailment, current/domain/power admission failure, finite SOC exhaustion, nonlinear refusal or cancellation publishes neither state. Motor power limits apply to the original held-discrete-voltage/linear-current reconstruction at the winding terminal, not total battery demand including converter loss and auxiliary demand or unresolved bus peaks. Original averaged electrical and constant-operating-point mechanical limits remain: no switching/protection events, thermal drift, stability/uniqueness or arbitrary-step accuracy claim. `rotordynamics_battery_drive` covers independent full physical equations and work, analytic mechanical regeneration into chemistry, and transactional refusal/retry.

## Coupled normal and torsional JKR motion on native bodies

`tangential_friction::adhesive::jkr_torsional::PreparedRigidJkrTorsionalPair`
connects the original radius-relaxed mode-I/III Griffith potential to actual
paired rigid bodies. Its two coordinates are normal approach in metres and
retained elastic twist in radians. The complete normal-point-force/free-couple
Gram matrix includes the off-diagonal response from actual lever arms and full
inertia tensors; no independent effective masses or changing-radius coordinate
substitution replaces it. Normal force, torsional moment and attached area
therefore evolve together through the original work-consistent source path.

Collinear original material markers and an immutable normal datum authenticate
approach. Reactions act at one common midpoint plane so internal force and
orbital-plus-spin angular impulse cancel. Normal-axis twist is source memory
driven by actual relative angular velocity, not an arbitrary difference of
endpoint quaternion logarithms. Held external COM forces and body torques
enter the same coupled solve. Contact and applied impulses are accumulated
before each body is updated once; compensated kinetic increments and actual
midpoint COM/angular travel audit independent external and contact work.

The native checks reconcile both mixed-unit coordinate increments and rates,
the normal marker/reaction projection, individual impulses, total linear and
angular momentum, and kinetic-plus-original-storage change against applied
work. Elastic release can assist body motion without producing fictitious
heat. Relative tangential point motion and transverse relative spin are
screened at both endpoints, which bound these affine rates throughout the
held-effort phase. The original source admits the entire quadratic trajectory
against its attached branch and material limits, including hidden reversals.

Preparation and every failed source, cancellation, nonlinear, geometry or
native audit leave both bodies and elastic history unpublished. A state pins
source chronology, native bodies and properties, body indices, original frame,
marker arms and normal datum. After normal travel, another held-pose phase
requires a compatible geometry continuation; the normal datum is never moved
to hide stale geometry. Sealed COM travel and body angular travel retain a
sticky eligibility flag for commuting principal-axis drift. This flag admits
only roundoff-level departures of each individual world angular velocity and
world spin momentum from the original normal, independent of user audit
tolerances. Cancelling a transverse excursion cannot erase it.

`continue_axial_geometry` consumes that sealed travel into actual COM positions
and exponential-map orientations, restricted to commuting spin about the
original fixed normal. It independently checks world spin, total orbital-plus-
spin angular momentum, full native kinetic storage, actual marker approach,
and both projected rates. It rebinds the original motion owner to the actual
new mobility while retaining the exact source coordinates, cached response,
time and step count. No caller attitude reset, coordinate projection, impulse,
extra time or compensating work is introduced. The original body momenta stay
unchanged, and the continuation returns both new captured pair and retained
state together. Noncommuting paths or unresolved orbital drift refuse.

This is reversible, small-deformation, ideal full-stick I/III motion on a local
attached branch. General noncommuting attitude paths, shear/rolling compliance,
partial slip, calibrated mixed-mode toughness, radial snap, detached landing,
reattachment and heat transport remain outside this owner. The focused
`adhesive_jkr_torsional` tests exercise changing area and original source work,
the zero-twist normal limit, an independently applied full-inertia impulse
oracle, off-axis held preload, actual applied work, elastic energy return,
world-rotation covariance, stale geometry and atomic unrepresented-mode and
cancellation refusal with retry. Geometry regressions exercise repeated actual
pose evolution with retained twist and clock, energy conservation, cancellation
and retry, and refusal of transverse finite rotation despite loose audit gates.

### Source temperature continuation of native normal/torsional adhesion

`PreparedRigidJkrTorsionalPair::continue_at_temperature` authenticates the
original held-indentation thermal refresh against the retained I/III source
and actual native history. The same original source re-equilibrates the disk at
unchanged approach and elastic twist. Signed parameter work includes normal
elastic, torsional elastic and surface storage changes at that relaxed area;
the refresh's pure-normal work cannot replace the twisted-contact work.

The immutable proposal retains both body poses and momenta, mass properties,
full mixed-unit mobility, original frame, markers, normal datum, coordinates,
rates, clock, steps, budgets and all pending geometry travel exactly. Only the
authenticated source and cached constitutive endpoint change. Native geometry
and source chronology must already match: normal travel requires compatible
axial geometry before a stale held-pose state can accept a temperature update.
An accepted pending twist path is retained through the update and remains
available to the original axial geometry owner afterward.

The source and native signed-work gates both apply before returning the new
pair and retained token. Source or applicability mismatch, invalid successor,
stale native data and cancellation publish no state. Temperature can supply or
remove stored mechanical energy without an impulse or fictitious friction heat.
This endpoint continuation supplies no caloric equation, thermal expansion,
continuous thermal-path stability, time increment, slip or reattachment.
Focused native regressions cover unequal bulk/interface temperatures, warming
and cooling work, changed subsequent recoil, no-op chronology, actual geometry
handoff, retained pending attitude drift and atomic refusal with exact retry.

## Passive full spatial-wrench relaxation between finite native bodies

`spatial_loads::internal_impulse::viscous::PreparedViscousConnector`
prepares the original full six-coordinate passive damping matrix at two fixed
native poses, material attachments, full inertias and one interval. The relative
twist is evaluated at attachment B in A's material axes, transporting A's rigid
velocity to B before subtraction. Common rigid translation and rotation are
therefore physical null motions even at separated attachment points. The
original internal-wrench owner transports the opposite free couple, so
transverse shear cannot supply an unaccounted orbital torque.

The complete relative mobility includes both finite masses, all off-center
rotational recoil and translation/rotation cross terms. Original physical
coordinate scaling and native Jacobi modes integrate the constant-geometry
passive linear source exponentially. Only the original passive matrix's nonzero
Gram-factor columns enter the at-most-six-dimensional eigenproblem; impulses
remain in that physical damping range. Stiff intervals cannot turn an exactly
absent damping direction into a numerically invented decay mode. No damping
floor, rank threshold or matrix projection is introduced.
Preparation uses one bounded native eigenproblem; retained state is fixed-size and
each relaxation allocates nothing, factors nothing and evaluates no exponential.
Actual pose changes require re-preparation.

One simultaneous original internal impulse updates both bodies. Original
per-body and total linear/angular momentum and work gates remain active.
An independent modal integral supplies nonnegative physical dissipation;
accepted kinetic change, connector work and recovered endpoint relative twist
must agree with that source to explicit dimensional tolerances. Numerical
defects are reported separately and never become heat. Source errors, stale
poses, numerical refusal and cancellation publish neither endpoint.

This is a composable frozen-geometry source step. It does not move poses,
impose a joint, add spring storage or allocate heat to a thermal reservoir.
Changing-geometry trajectories require composition with the original native
motion and temporal refinement. The material thermal axial dashpot and the
rotation-only viscous pair retain their existing source and lifecycle owners.

Focused tests in `tests/internal_viscous_connector.rs` cover analytic weak/stiff
axial recoil and heat, transverse shear with its necessary transported couple,
all six coupled modes against independent continuous Newton/Euler integration,
full inertia, physical coordinate rescaling, common-rigid-motion nullity,
cross-coupled rank-one free modes under extremely stiff exposure,
distant origins and common boosts, zero damping, stale geometry, cancellation
and undeliverable native momentum transfers.

### Simultaneous held actuation and full spatial damping

`viscous::driven::PreparedDrivenViscousConnector` integrates a held internal
force/free couple and the original full passive damping on the same fixed-pose
trajectory. With native relative mobility `G = B^T M^-1 B`, the conjugate
coordinate `x = G^-1 s` obeys `x' = f - C G x`. The existing pHS exact-ZOH
owner propagates all six coordinates, integrates signed actuator work and
independent positive heat, and retains the original loss-factor witness.
This supplies the physical terminal slip without alternating force kicks and
unforced damping. Semidefinite and zero damping retain driven free directions.

One native transported internal impulse applies the accepted conjugate momentum
change to both original finite bodies. Initial coordinate reconstruction,
actual endpoint twist, original per-body and total momentum, native work and
independent heat/work balance all gate publication. Preparation is bounded to
six relative coordinates; warm calls allocate nothing and reuse their original
factor and exponential/heat caches. Poses, inertia, damping and interval are
bound; input momenta and the held effort may change between calls.
This source does not advance poses, introduce stiffness or assign a thermal
destination. Changing geometry still needs temporal refinement. G1 tests in
`tests/internal_viscous_connector/driven.rs` check scalar terminal velocity,
actuator energy recovery, independently integrated heat, full spatial recoil
against continuous Newton–Euler integration, subdivision, driven null losses,
zero/unforced limits, stale geometry, cancellation and deterministic reuse.

## Finite motion with passive full spatial damping

`spatial_loads::internal_motion::viscous` advances two original finite bodies
through exact fixed-geometry viscous half sources around one native loaded
rigid-body step. The endpoint source is prepared from the actual endpoint
poses, material axes and attachment transport, with the same authored
`PassiveBushingMatrix`. All six coupled force/couple coordinates and the
original complete inertia tensors participate. The native free-motion,
held-force/body-couple, and prepared material-point/world-couple paths remain
the respective pose and external-work owners.

`step_viscous_connector`, `step_viscous_connector_with_loads` and
`step_viscous_connector_with_spatial_loads` return actual changed poses and
momenta, both source receipts, independent viscous dissipation, external
force/couple work, gravity work/potential change, and whole-step momentum and
mechanical-energy/heat residuals. The original retained kinetic-storage
reduction supplies actual quadratic changes; heat is never inferred by
subtracting observed energies. Each native drift and each viscous source
passes its own balance gate before the complete original-to-final account.
Numerical drift defects remain separate from physical dissipation.

The smooth frozen-coefficient composition is second order under time
refinement. Exact local source stability does not establish changing-geometry
accuracy or remove the need to refine stiff trajectories. Work consists of
two bounded six-coordinate preparations and two native body steps, including
the original bounded spectral scratch and native nonlinear/rotation limits.
No hidden subdivision, new inertia owner, projection, growing state history,
joint/contact law or thermal allocation is introduced. Inputs are immutable;
invalid duration, source refusal, native solve/balance refusal or cancellation
returns no partially accepted endpoint.

Focused tests cover analytic finite-mass recoil and heat with a three-rung
position refinement, a full moving six-axis source against independent
continuous Newton-Euler RK4 dynamics, native zero-law equivalence, common
boost/world-translation covariance, original rotating material-point work,
gravity/heat separation, and late cancellation with exact retry.

## Finite motion with simultaneous internal actuation and spatial damping

`spatial_loads::internal_motion::viscous::driven` advances the original pair
through simultaneous driven/damped half sources around one native loaded
rigid-body step. The complete passive matrix and full finite-body mobility act
on the same source trajectory as the actuator. The endpoint source is prepared
from the actual evolved poses, material axes and attachment transport.
`HeldInternalEffort` explicitly distinguishes world-held force/free-couple
components from components following body A's material axes; world components
are re-expressed at each actual source geometry.

Each source independently integrates signed actuator work and nonnegative
viscous heat. Native COM/body-couple, gravity and rotating material-point loads
retain their original pose, work and impulse owners. Every source and native
body drift passes its own gate before the original-to-final momentum and
mechanical-energy account. Numerical drift defects remain separate from heat;
actuator energy recovery is signed work. Cancellation or refusal publishes no
partial endpoint, and no thermal destination is inferred.

Two bounded six-coordinate preparations are released in sequence around the
native body steps. This smooth composition is second order under temporal
refinement; local exact stiff-source propagation does not certify coarse
moving-geometry accuracy. G1/G3 regressions in
`tests/internal_viscous_motion/driven.rs` cover analytic terminal slip,
independent work and heat, refining axial positions, both actuator frames against
continuous full-spatial Newton-Euler integration, passive/zero-law limits,
native external-load work, gravity separation and cancellation with exact retry.

## Sparse full-spatial viscous connector networks

`spatial_loads::internal_impulse::viscous::network::PreparedViscousNetwork`
binds original poses, full mass/inertia, attachment arms, passive six-coordinate
laws, source duration and an explicit positive subdivision count. Arbitrary
shared-body graphs, cycles, parallel connectors and isolated bodies are supported.
Each active edge reuses the original exact finite-body connector, including
material-frame force/couple cross-coupling and opposite-wrench moment transport.

A forward-half/reverse-half composition merges the central active edge into one
full exchange. Each subdivision uses at most twice the active edge count minus
one native pair relaxation. Zero-duration and zero-law edges perform no pair
work; per-edge heat remains indexed in the original supplied order. The retained
responses perform no source-loop eigensolve or exponential evaluation. Preparation
admits body, connector, pair-call and logical peak storage limits before allocation,
including the original bounded small-eigenproblem scratch. Source work scales with
active edges times subdivisions, with linear body/edge admission and readback.

The composition is second order under subdivision for this fixed-geometry source;
edge order is part of the finite-step method. Local stiffness stability does not
certify a coarse simultaneous-network solution. No pose, clock, constraint,
elastic storage or thermal state is advanced. Reprepare after any original
geometry, inertia, law or interval changes.

All source updates remain in a private native state bank. The original pair
receipts admit each local physical response; a final all-body audit additionally
compares accumulated unrounded intended force/moment transport with actual
represented body increments. Total linear and orbital-plus-spin angular
momentum close about the caller's declared origin. The original
`RigidKineticEnergyChange` retains all bodies' endpoint products before
rounding; both that kinetic increment and independently accumulated connector
work must balance analytically integrated nonnegative dissipation. Per-connector
heat is reported for explicit thermal coupling. Refusal or cancellation exposes
no intermediate accepted body states.

The underlying pair preparation exposes
`maximum_decay_rate_per_s()` for downstream material applicability
and subdivision admission. It retains the largest original squared kinetic mode
rate without another eigensolve. Zero-duration or zero-law preparation returns
zero because no active spectrum was prepared; this is not a claim that a
nonzero constitutive law has zero physical decay.

Focused regressions cover exact single-edge equivalence, analytic three-body
graph modes, shared cycles with full nonprincipal inertia against independently
integrated force/torque/heat equations, quadratic refinement, common rigid and
isolated free modes, distant origins and common velocity boosts, inactive
sources, explicit resource boundaries, late physical refusal, and cancellation
with unchanged-input retry.

## Moving full-spatial viscous connector networks

`spatial_loads::internal_motion::viscous::network` composes the original passive
graph sources with one native finite-motion step per body. Shared bodies,
cycles, parallel connectors and isolated bodies retain their complete mass and
inertia. Both half sources are prepared from their actual endpoint geometry;
evolved attitudes, off-center arms and transported opposite couples therefore
participate in the second exchange. Held native loads, uniform gravity and
prepared material-point/world-couple loads use the existing motion/work owners.

`ViscousNetworkMotionControls` separately authors source subdivisions and
physical balance tolerances. Its resource limits cover both stages' total pair
calls and the complete logical peak storage before allocation. Each source
preparation is released before the next, and the final source receipt owns the
endpoint state bank without a duplicate output copy. Storage stays linear in
bodies plus connectors; no global matrix or new source kernel is introduced.

The returned per-connector heat is the sum of independent positive source
integrals, available for explicit thermal allocation. Applied work, gravity
potential, native numerical drift defects and physical damping loss remain
distinct. Every native body drift passes its own gate before the complete
original-to-final momentum and mechanical-energy/heat audit. Any source,
motion, resource or cancellation refusal publishes no partially accepted state.

The smooth moving-geometry composition is second order under temporal
refinement. The graph source is separately second order under its explicit
subdivisions; passive exact pair exchanges do not make a coarse whole-graph
trajectory exact. No joint constraint, spring storage, temperature evolution,
thermal coefficient update or microscopic damping mechanism is inferred.
`tests/internal_viscous_network_motion.rs` exercises G1 continuous graph-mode
position/velocity/heat refinement and a full moving six-axis cycle against
independent Newton-Euler RK4; G3 original pair/zero-law/isolated-body behavior,
gravity separation and rotating point-load work; and bounded-resource and late
cancellation refusal with exact retry.

## Native simultaneous full-stick normal, shear and twist adhesion

`adhesive::jkr_combined` consumes the original four-coordinate circular JKR
source on one actual pair of finite bodies. Its coordinates are normal approach,
two elastic tangential displacements and elastic twist in radians. One relaxed
contact area supplies all four conjugate efforts and normal, shear, torsional
and surface storage; independent normal/shear/twist disks, Coulomb caps and
radius-scaled twist histories are not substituted.

The complete mixed-unit `J M^-1 J^T` comes from actual masses, full inertia
tensors, the original frame and coincident reaction arms. Source efforts act as
opposite point forces and opposite free normal couples. External held COM forces
and body torques enter the same coupled solve through their actual accelerations.
Existing native impulse owners apply the combined resultant once to each body;
independent point/angular rates, integrated travel, actual kinetic change,
external work, source storage and body/pair linear/angular momentum must agree
before both momenta and sealed constitutive history can be published atomically.
Returning reversible contact energy is assisting mechanical work, never heat.

Only elastic shear/twist preload is declared at initialization: material-marker
geometry supplies approach and actual body momenta supply all four rates.
History retains original source, mobility and numerical-policy identity together
with body indices, mass properties, poses, momenta, frame and material datum.
The phase holds geometry and retains accumulated COM and body-rate travel; those
integrals are not finite attitudes. A changed normal approach refuses subsequent
use against stale marker geometry. This owner supplies no generalized pose
continuation, slip, rolling elasticity, radial snap or detached trajectory.
Relative transverse rolling is bounded at both endpoints and therefore throughout
the affine held-resultant phase. Common transverse spin is allowed and contributes
its actual contact-point velocity; facing equal-inertia bodies can consequently
recoil under both shear directions without discarding their point-force torques.

Five focused regressions exercise simultaneous four-coordinate recoil and external
work with changing area, all sixteen full-inertia mobility entries against actual
unit body impulses, off-axis loaded equilibrium, elastic energy return and world
rotation covariance, original torsional/normal limits, and atomic original-source,
body, stale-geometry, rolling and cancellation refusal with deterministic retry.

### Physical rotor modes, residues and shaft-speed slopes

`RotorOperatingPoint::prepare_resonance_modes` retains the complete original
quadratic pencil `s² M + s(C + Omega G) + K`. The existing `fs-modal` owner
supplies binary coordinate/rate equilibration and local simple-pole refinement,
with original-coordinate unit-mass displacement and velocity shapes and source
residuals. Signed decay, oscillation and active growth remain intact. Existing
all-pole and Floquet queries keep their separate scope.

`PreparedRotorModes::prepare_response` rechecks the mode against the current
original pencil and retains one shared bordered factor for true nonsymmetric
force/receiver residues and physical parameter derivatives.
`PreparedRotorPole::spin_derivative` supplies the analytic direction
`dT/dOmega = s G`, holding M, C, K and per-speed G fixed; it also differentiates
the normalized displacement and velocity shapes. Speed-dependent bearing,
preload or material coefficients require their own actual law.

The rotor convention is `exp(+i omega t)`: complex angular frequency is
`-i s` and its Laurent residue is `-i` times the Laplace residue. Force
participation uses the full bordered equations, not a right-mode projection.
Zero forces, coherent cancellation and dark receivers can return zero.
Physical modes and one-pole residues do not substitute for a driven peak,
complete modal expansion, crossing identity or stability certificate.
Repeated/defective roots, resource limits and cancellation retain the shared
owners' explicit refusals.

### Hereditary rotor-bearing resonance

`RotorBearingModel::at_speed_with_maxwell_memory` retains original scalar strain
channels, positive branch stiffnesses and relaxation times at a signed operating
point. The supplied K is the relaxed stiffness. Each original branch obeys
`sigma_dot + sigma/tau = k*b^T*q_dot`; its unrelaxed stiffness must not already
be added to K. Original nonsymmetric/active C and K and signed gyroscopic forces
remain in the mechanical equations.

Harmonic response uses the exact original Maxwell transfer and factors only
the mechanical n coordinates, at O(n^3+b*n^2) work. The same LU supplies complete
source/frequency/spin derivatives. Every original branch stress, cycle-mean
stress work, passive heat and constitutive residual is recovered, and mechanical
equilibrium plus signed power balance are separately checked. A vanished
nonzero-frequency loss quadrature refuses rather than silently dropping heat.
No sampled loss factor or damping-ratio fit replaces the branch law.

`poles` delegates the complete 2n+b state spectrum to the shared Maxwell owner,
retaining equal-rate dark and uncoupled material poles. `prepare_modes` retains
the shared local mechanical-mode/residue domain, which excludes pure material
poles and the relaxation denominators. `driven_resonance` resolves the selected
physical response peak and optional half-height width using original dispersion.
Branch count is capped at64; native cancellation and resource refusals remain
explicit. Coefficients describe the supplied operating state, without inferred
preload, temperature or speed-dependent constitutive laws.


## Hereditary synchronous and unbalance resonances

`RotorMaxwellOperatingPoint::synchronous_response_with_derivative` follows
`omega = order * Omega` for finite nonzero signed order and signed shaft speed.
It retains source phase, the separate gyroscopic direction, and every original
relaxation derivative in the same native factor solve. `synchronous_resonance`
samples the same immutable model and branches over positive shaft speeds. Peak,
brackets, half-height width and observation derivative use **shaft-speed** units,
including when order differs from one.

`RigidRotor::unbalance_maxwell_response_with_derivative` and
`unbalance_maxwell_resonance` reuse original station, phase and moment forcing
`F = U * Omega^2` and its total derivative, with zero source and derivative at
rest. Original branch stresses, passive heat and physical residuals remain in the
returned response. No frozen-force or frozen-loss approximation is introduced.
The existing model stiffness is the relaxed stiffness; branch stiffness is not
counted twice. These compare steady operating points with supplied coefficients,
not startup history or inferred bearing/preload/temperature laws.

`rotor_maxwell_synchronous` supplies five independent native test groups,
including a rational hereditary unbalance peak and both half-height roots,
signed-order shaft-unit scaling, source moments, cancellation and refusals.

## Native rolling through adhesive separation and recontact

`tangential_friction::adhesive::rolling::events` composes one original
prescribed-normal lifecycle stage with an unforced native midpoint rolling
phase. Original JKR attachment/separation thresholds remain authoritative.
Attached endpoints use both actual inertia tensors and the original source
moment solve; its accepted angle is replayed once through the sealed lifecycle
owner to retain cumulative release and plastic-loss accounts. A detached pair
keeps its original momenta and records bounded free angular travel without
preloading a future junction. A normal snap supplies no invented body impulse.

Preparation binds the source history to the original pair, material arms,
poses, mass properties and frame. It compares original mechanical endpoints
with the regenerated lifecycle, checks the separate incremental energy channels,
and independently observes both normal and rolling endpoint storage with an
explicit rounding allowance. All fallible checks precede simultaneous
publication of the two momenta and full lifecycle state. New native bindings
create fresh attached or detached histories; loaded memory cannot be silently
adopted into another frame. G1/G3/G4 tests in `tests/adhesive_rolling/events.rs`
specify analytical spring release, absence of a spurious snap impulse, exact
fresh-contact recoil after reattachment and atomic late-refusal/retry behavior.
The scope is prescribed indentation followed by unforced rolling at held poses
and frame; normal body motion, reciprocal area relaxation, heat destinations
and finite rotation remain with their separate physical owners.

## Native bonded torsion through prescribed normal loading

`tangential_friction::adhesive::torsion::normal` changes the original attached JKR radius at fixed retained twist, then runs the original native midpoint torsion solve with actual full body inertias and held COM forces/body-frame torques. Normal JKR work, signed torsional-stiffness work, applied body work, recoverable twist and calibrated plastic loss remain distinct. Original geometry, material/calibration identity, mode-III applicability, hidden-turn, momentum, bounded-effort and mechanical work gates remain active. Original incremental work primitives drive the combined balance; independent absolute-store observations receive only an explicit roundoff allowance. Opposite stage defects cannot cancel.

`prepare_rigid_adhesive_torsion_with_normal_update` stages the original source, retained native state and both bodies. The corresponding apply operation publishes all owners only after every fallible check and the final cancellation poll. This is prescribed normal loading at held poses/frame with homogeneous bonded twist; it supplies no normal impulse, detachment, stress-free annular recruitment, reciprocal Griffith equilibrium or thermal destination.

Three added G1/G3/G4 integration groups cover analytic loaded-spring recoil at tensile/compressive targets with actual anisotropic inertia, independent native kinetic/impulse checks, exact replay of unchanged-area wrench-driven motion, and atomic cancellation/source/geometry/wrench refusal with deterministic retry. Formatting and recursive syntax/module resolution passed; these integration tests were not type-checked or executed because DSR/RCH were unavailable and the existing full dependency build had exhausted shared memory/disk resources. The seven underlying production torsion leaf tests passed.

## Applied loads across adhesive rolling attachment and separation

The rolling lifecycle's `prepare_rigid_adhesive_rolling_events_with_wrenches` and matching apply entry point accept actual held COM forces and body-frame torques in original A/B order. Attached motion uses the original midpoint moment solver, including calibrated plastic loss and exact source-history replay. Detached motion applies each original resultant-impulse owner once, retaining actual represented momenta, compensated kinetic change, midpoint load work and per-body work/momentum defects. Summed world linear and orbital-plus-spin momentum residuals are checked independently. No normal snap impulse or detached rolling preload is invented.

Reports expose external work and actual kinetic change for both branches, with original native impulse receipts for separated bodies. Prescribed normal/rolling-parameter work remains distinct from applied-body work. Detached constant-acceleration turning excursions are checked against the angular bound even when the net angle is zero. The original zero-load entry points delegate to the same implementation. Both bodies and source lifecycle publish only after all gates and the final cancellation poll.

Three new G1/G3/G4 groups cover identical original driven attached behavior with plastic yield, independent detached force/torque momentum/work/travel and fresh recontact, zero-wrapper replay, invalid/overflow load rollback, and hidden-turn refusal. Pinned formatting and recursive syntax/module resolution passed; these integration groups were not compiled or executed under the documented build-resource limits. Poses/frame and prescribed normal motion retain the existing approximation; finite orientation, automatic normal-body motion and thermal conversion are separate owners.


### Native equilibrium with constant joint locks and gears

`articulated::motion::spatial_loads::prescribed::equilibrium::constrained` solves the original supported-tree force equations and constant `JointConstraint` rows simultaneously: `R(q)-u-A^T lambda=0` and `Aq=b`. Stationary and symmetry-admitted uniform-rotation paths use the same native RNEA/world-load/finite-attachment residual as motion. Smooth stationary actuator feedback retains the original repeatable-trial contract. Ideal reactions remain separate from actuator limits, storage, dissipation and support work; no penalty spring or relaxation trajectory is invented.

Explicit joint, effort, multiplier and row-position scales define a bounded matrix-free Newton-Krylov problem through the same dispatch as unconstrained equilibrium. The original at-most-32-row mobility owner checks dimensions and normalized independence at the seed and accepted root. Every trial retains original joint/source admission. Accepted force residuals use the actually representable native row reaction arithmetic, and native absolute position tolerances are independently checked after solver convergence. All residual work, final admission and cancellation remain bounded; refusal leaves the seed unchanged. Empty rows delegate exactly to the original unconstrained solver.

The result carries actual actuator and ideal-reaction efforts, multipliers, force and position residuals, original support wrenches and attachment storage. It initializes the existing constrained supported-motion owner. Contact, backlash, changing ratios and nonlinear geometric closure are outside this affine model. Stability, uniqueness and constrained small-signal response are not inferred from force balance; the unconstrained linearization cannot replace a tangent that retains these same rows. Regression cases cover loaded geared springs and native motion continuation, a fully locked zero-stiffness joint with limited actuator capacity, signed row rescaling, uniform spin, actuator feedback, empty-row parity and bounded refusals.

## Native bonded torsion through separation and fresh recontact

`tangential_friction::adhesive::torsion::events` advances the original sealed torsion lifecycle and both finite rigid-body momenta together. The source's prescribed monotone normal stage owns attachment, the actual mode-I fold, normal work and normal snap release. Retained twist follows the changing area to that fold; its remaining spring energy exits separately from plastic torsional loss and supplies no body impulse. Fresh recontact starts an unstrained spring while preserving cumulative loss and release histories.

Attached motion uses the original shared finite-inertia angular solve and checks its exact constitutive replay. Actual held COM forces and body-frame torques act on both attached and separated pairs, with original momentum receipts, midpoint work, native kinetic-energy changes and bounded angular travel. Separated travel cannot preload a future junction. Rolling and torsion now share the private detached-body phase without changing its arithmetic. All source and native endpoints are staged before the final cancellation check and publish together.

The normal, stiffness, spring-release, plastic-loss and body-work channels remain distinct, with compensated incremental accounting and independent stage gates. Poses and frames stay held; this retains the original homogeneous twist law and explicit mode-III/adhesion applicability ceiling. It does not supply reciprocal mixed-mode area stability, annular recruitment, normal body integration or a caloric destination.

Four focused integration cases cover the independently derived fold-storage ratio, exact driven/changing-area native replay, detached force/torque work and fresh recontact, and seven atomic refusals followed by deterministic retry. Pinned `rustfmt --check` passed. These integration cases were not compiled or executed in this constrained session; no new runtime-verification claim is made.

## Simultaneous driven spatial damping on finite graphs

The original viscous-network module now exposes a separate
`network::driven::PreparedDrivenViscousNetwork` source. Every native body
retains its six kinetic coordinates and every connector retains its original
material-frame full passive matrix. The original transported native wrench
columns produce the work-conjugate map B. The complete system
`y' = B f - B C B^T y`, `H = y.y/2` goes through the existing exact pHS
ZOH owner; cycles, parallel links, isolated bodies and undamped directions
need no inverse edge mobility, rank cutoff or ordered pair approximation.
The existing sparse palindromic passive source is unchanged.

Each original connector has an independently integrated positive source-square
heat observer on that same complete trajectory. Signed held-actuator work is
reported per connector and in total. Supply and recovery remain work; neither
numerical drift nor an energy-account remainder becomes heat. Per-link heat
and work sums must agree with the full pHS integrals. Integrated net efforts
enter the existing simultaneous native impulse network once, preserving original
force/couple transport, full inertias and actual native P/L/work admission.
Initial kinetic-coordinate reconstruction and actual endpoint relative twists
are checked against the source trajectory before publication.

The moving graph entry reuses the original loaded native body motion between
two complete graph half sources. Material-held effort follows the first body;
world-held effort is re-expressed at actual endpoint geometry. Applied loads,
gravity storage, signed actuator work, material heat and native drift defects
remain separate. Every source, each body drift and the complete original-to-final
account must pass. Borrowed input states, laws and efforts never change on
refusal or cancellation.

Explicit body, connector, peak-byte and preparation-work caps cover all retained
source observers and both moving preparations before allocation. Both the doubled
body dimension and the state-plus-input integral dimension must fit the original
256-coordinate pHS bound. Warm source calls use cached contractions, with no
exponential or factorization; the original native event and returned body/link
banks allocate bounded linear storage. Moving geometry still needs independent
temporal refinement. No contact constraint, spring model, thermal destination or
generic DrivenBody description binding is inferred.

Six fixed-source regression groups passed natively in an isolated build of the
actual source and original numerical dependencies. They check analytic cyclic modes and
terminal slip, simultaneous positive/recovered work on parallel links, full
tensor pair parity, an independent native-coordinate force/torque/heat ODE,
weak loss beside a stronger source, and atomic source/budget refusals. Three
moving groups specify original loaded-pair parity, second-order positions with
exact fixed-geometry graph rates, and complete-step refusal/retry. Pinned Rust
formatting passed. The focused source run is not full workspace integration;
full moving/thermal integration execution remains outstanding.


### Ideal lock/gear equilibrium response

`articulated::motion::spatial_loads::prescribed::equilibrium::constrained::linearized`
re-admits original affine positions/rates, native row mobility and RNEA balance
against actual A^T lambda. Cached equilibrium reports do not admit a root;
the unconstrained response retains its no-reaction gate. Original Householder
QR and shared fs-la subspace arithmetic form a complete physical lift dq=N z,
independently checking original A N and scaled orthogonality. Only allowed
directions enter native K/C stencils, with work-conjugate M/C/K projection.
Gyroscopic nonsymmetry and negative stiffness remain in the existing modal
owner. Fully locked/empty models return an empty pencil without multiplier
modes or derivative probes. Explicit source work, cancellation and range gates
refuse without changing original state. This path admits constant affine rows;
nonlinear closures, changing ratios and contact transitions need their own
complete tangents. Five analytical native integration regressions are authored
but unexecuted here; shared subspace arithmetic passed eight original-source
checks. Reconstructed integration bytes received source review after executor
loss and do not claim a completed full-crate build or stability certificate.

## Forced periodic rotors with complete bearing memory

`periodic::{prepare_maxwell_forced_cycle, prepare_maxwell_harmonic_cycle,
prepare_maxwell_fourier_cycle}` extend the same complete original
`[q, q_dot, sigma]` cycle with held, harmonic, or simultaneous DC/Fourier
generalized force/torque. Loads enter mechanical momentum rows only; every
original stress responds through the complete reciprocal generator. Full M/G,
signed spin, branch strain/stiffness/time and relaxed-K conventions are unchanged.
Source phase, displacement, velocity and all stresses remain continuous through
declared coefficient and load-amplitude switches. A load step is not an impulse.

The existing normalized affine/harmonic exponential and Fourier superposition
owners now accept the admitted full physical state dimension. Their original
mechanical callers still use 2n; hereditary callers use 2n+b. A harmonic source
adds two auxiliary coordinates (at most 130 total) and held forcing adds one;
none enter the original physical Floquet spectrum. Fourier offsets accumulate
in double-double arithmetic while retaining one homogeneous period map. One
mass factor and one original augmented generator per segment serve all tones.
Budgets count every declared source and full state before dense assembly;
retained maps and warm advance/orbit costs are independent of tone count.

`PreparedRotorMaxwellCycle::forcing_offset` exposes the original physical
affine offset and `periodic_orbit` returns `RotorMaxwellPeriodicOrbit` with
every boundary stress, original fixed-point residual, independently advanced
closure residual and balanced resolvent-gain estimate. Mechanical and hereditary
orbits share the existing native LU/condition-admission owner. Singular families
and unresolved slow directions refuse even when their loads vanish; no state is
removed or phase pinned. A unique unstable particular orbit retains its actual
positive Floquet growth and does not imply attraction.

Eleven focused regression groups are authored for exact three-exponential held
response and asymmetric load cycles, relaxed static equilibrium, dark/slow
memories, full inertia and signed gyroscopic congruence, original complex
Maxwell phases/stresses, DC plus tones, independent original-ODE transients,
secular lossless forcing and cancellation/admission. Independent source and
physical-oracle review was followed by native execution: all eleven groups now
pass in the focused actual-source Cargo build, together with six existing
homogeneous-memory groups and the ordinary forcing regressions. The build uses
the original numerical modules without substitute implementations; it is not
a full-workspace or DSR verdict.

The declared piecewise coefficients and finite source spectrum carry no smooth
modulation or waveform-truncation error bound. The accounted counterparts below
add physical work and heat to these same trajectories. Nonlinear bearing laws
and eigenvalue enclosures remain separate capabilities.

## Work, storage and heat of forced Maxwell rotor cycles

`periodic::{prepare_maxwell_forced_power_cycle,
prepare_maxwell_harmonic_power_cycle, prepare_maxwell_fourier_power_cycle}`
prepare the original full-state motion together with independently integrated
physical work channels. `PreparedRotorMaxwellPowerCycle::advance` handles
arbitrary initial displacement, velocity and branch memory; `periodic_orbit`
accounts the native orbit's already-computed trajectory without a second
propagation. Original fixed-point, closure, resolvent and Floquet diagnostics
remain available through the motion owner.

Each `RotorMaxwellSegmentWork` separates applied generalized-force work,
signed `integral v^T C v dt`, signed circulatory `integral v^T K_skew q dt`,
kinetic storage, relaxed symmetric-K storage, every original
`Delta(sigma_j^2/(2 k_j))`, and each independently integrated branch heat
`integral sigma_j^2/(k_j tau_j) dt`. Dark branches and equal relaxation rates
retain their original identities. Gyroscopic force remains in motion and
contributes no heat; negative active-bearing and circulatory work remains
signed work. A negative numerical Maxwell heat refuses rather than clipping.

C and relaxed K may change at declared boundaries with continuous q/v/sigma.
Outgoing stiffness work is `q^T(K_next-K_old)q/2`, including the final-to-first
period seam. Cycle endpoint storage uses the first segment's coefficients on
both sides. The complete balance is storage change plus C/circulatory work
plus branch heat, minus applied work and stiffness-switch work. Each segment,
the sum of absolute segment defects, and the complete cycle must independently
satisfy the caller's absolute energy allowance. `RotorMaxwellCyclePower`
retains all channels, per-branch heats, storage changes and mean powers.

Native motion and work assembly share the extracted original
`MaxwellGenerator`. The existing pHS quadratic-flow integral owner observes
the complete joint physical/source state; the existing source-phase congruence
then condenses source coordinates into one affine slot. DC/tone, tone/tone and
initial-memory interference remains present. Separate steady-tone powers are
not substituted for the switched transient. Original mechanical storage,
source insertion, switch-work and energy gates are reused without a separate
trajectory integrator or material-memory reduction.

The original fixed M/G, signed spin, branch strain vectors, stiffnesses and
relaxation times remain required. Changing these branch storage laws needs
separately accounted parameter work. The returned heat has no inferred thermal
destination or temperature law. Closure is a numerical energy account, not an
interval trajectory enclosure, passivity proof or orbit-attraction guarantee.
Smooth coefficient modulation and finite Fourier truncation still require
their own convergence evidence.

With s segments, b branches and h tones, cold work is
O(s*(3+b)*(2n+b+1+2h)^3); retained forms and warm accounting are
O(s*(3+b)*(2n+b+1)^2). Complete native motion, concurrent quadratic workspace,
retained forms and result storage are admitted before preparation. At the
original hard limits the joint generator has at most 193 coordinates, within
the shared 256-coordinate quadratic owner. Factors/forms are immutable after
preparation and cancellation publishes no partial account.

All six new physical groups pass natively: independent dark-memory and coupled
decay heats; resonant periodic phasor consumption at source scales
1e-60, 1, 1e60; full-mass/gyro/circulatory coordinate invariance; switched Fourier
transients against an independent original-coordinate RK4 work/heat reference
with seam work and harmonic cross terms; zero-branch equivalence with the
existing mechanical power owner; and admission/cancellation/replay.
The focused actual-source validation totals 71 passes including the motion,
lossless-Fourier and notch/peak sensitivity regressions. DSR and RCH are
unavailable in this environment; these local Cargo runs are a bounded native
fallback, not a full-workspace or DSR gate.

## Open physical-time Maxwell rotor intervals

`periodic::{prepare_maxwell_interval, prepare_maxwell_power_interval}` advance
an authored open window at absolute time. `RotorIntervalForcing` supports
homogeneous, held, harmonic and finite multi-tone forcing with explicit
positive angular frequencies [rad/s], independent of the window length or
shaft spin. Noncommensurate tones are admitted in strictly increasing order.
Global source phase survives coefficient, amplitude, clock and window cuts;
changing a thermal timestep does not change the excitation frequency.

The existing Maxwell generator, source insertion, affine propagation and
quadratic-work owners are shared through one retained flow. Cycles alone add
their genuine QR spectrum and periodic-orbit operation. The distinct interval
types expose advancement and original physical transition/residual reports,
without a Floquet or attraction claim. Interval admission skips QR sweep
requirements and iteration work while retaining mass/assembly, exponential,
composition, trajectory and work-form bounds. Strongly decayed directions do
not require a representable Floquet multiplier merely to advance a state.

Motion and heat use the same DD absolute start-plus-elapsed source clock.
Physical phase uses the original binary frequency/time inputs and split-2*pi
reduction, with |omega*t| <= 2^40 radians and at most 1024 rotations per segment.
Rounded public boundaries must be distinct; their rounding relative to the
authored adjacent duration, and the induced phase error in radians, must fit
the caller's generator tolerance. Underresolved time handoffs refuse rather
than silently resetting source phase. Existing cycle phase arithmetic is
unchanged. The allocation-free forcing `validate` reuses the original bounded
shape admission and also supports composing owners' zero-duration checks.

The shared trajectory has absolute boundary times for intervals. Its power
account uses the actual first/last storage charts and only INTERNAL relaxed-K
switch work, with no final-to-first seam. If a caller changes K between separate
windows, that external parameter switch needs its own work account. Fixed M/G,
signed spin, original branch stiffness/strain and complete stress identity
remain required; optional held kinetics derive from original reference times.
Every original branch heat, signed C/circulatory work and force interference
remains present. Smooth coefficients and changing-temperature feedback still
require their own time refinement; numerical residuals are not error enclosures.

Four new native groups pass: motion, work and branch heat against independent original-coordinate
RK4 integration under switched clocks and noncommensurate tones;
open endpoint storage and internal switch work; unequal-window phase/work
equivalence and the 90-digit original-binary omega=1.3, t0=1e9, h=0.5 free-motion
oracle; complete dark-memory decay without QR; and shape, clock-resolution,
phase, resource and cancellation refusals. All 41 affected original periodic
unit and harmonic/Fourier/Maxwell regressions also pass (45/45 total).
The prior strict unresolved-force and increment-corruption refusals remain.

These are focused native Cargo checks using the actual numerical sources.
DSR/RCH are unavailable; no full-workspace or DSR gate, trajectory enclosure
or continuously coupled thermal accuracy is claimed by this interval layer.

## Held Maxwell clocks throughout rotor cycles

`periodic::RotorMaxwellKinetics` supplies one positive relaxation-time ratio
for every original segment and branch. Its `floquet`, `prepare_cycle` and
`prepare_power_cycle` methods reuse the existing full-state preparation,
native exponential, affine orbit and quadratic-work owners. The public
`RotorPeriodicForcing` selects homogeneous, held, harmonic or finite Fourier
loads with the original common-cycle phase convention. Homogeneous power
accounting uses no fictitious source coordinates.

Every effective time is `tau_reference * ratio[segment][branch]`, derived
independently from the original branch. Previously shifted times are never
compounded. The same admitted effective time enters the segment generator,
independent original constitutive residual and branch heat form
`sigma_j^2/(k_j tau_j)`. The original mass, signed gyro/spin, branch stiffness,
strain direction and complete stress identity remain unchanged. Stress stays
continuous at clock switches and the final-to-first seam. Changing only this
dissipative rate changes no stored elastic energy and requires no parameter
work; authored relaxed-K switches retain their existing signed work.

The allocation-free `validate` method checks the complete bounded layout,
positive finite ratios and positive representable effective times AND rates,
with cancellation per segment and branch. Dense preparation retains its
existing whole-query budgets; there is no extra exponential, material-pole
elimination or integration path. The schedule is borrowed only during cold
preparation; warm maps and work forms are immutable and need no subsequent
kinetic evaluation. Existing named fixed-clock entries use the identical
owner without a schedule.

These are caller-declared held kinetics, not an inferred temperature law or
thermal feedback model. Smooth rate variation still needs independent
time/coefficient refinement. Changing branch stiffness, strain coordinates or
topology requires its own history-transfer and storage-work law. Numerical
Floquet, trajectory and energy diagnostics retain their original limits.

Four new native groups pass: exact per-branch dark-memory decay and heat
across three noncompounded clock rows; coupled periodic Fourier response and
each branch's heat against independent original-ODE RK4 integration; bitwise
identity with the original motion/power owners at all-one clocks; and malformed,
unrepresentable, budget, cancellation and replay cases. All 25 affected existing
periodic-unit, Maxwell harmonic and Maxwell power regressions also pass (29/29).
The positive-time/infinite-reciprocal refusal is exercised separately from time
underflow. These are focused actual-source local Cargo runs; DSR/RCH remain
unavailable and no full-workspace gate is claimed.

## Resolved lossless motion at cancelling cycle boundaries

The ordinary and hereditary affine advance checks retain the largest actual
boundary-to-boundary increment separately in each original q, velocity and
stress coordinate. It supplements the original net-increment, force-offset
and (P-I)x scales when complete-cycle cancellation makes the net increment
nearly zero. Strong activity in another coordinate and a large unchanged
state never enlarge this allowance. Fixed-point queries supply zero activity
and retain their strict original force scale. The native trajectory, damping,
source coefficients and caller tolerance are unchanged.

The lossless Fourier fixture with simultaneous DC, resonant and off-resonant
forcing now advances with its analytic secular growth at amplitudes
1e-80, 1 and 1e80. Both new activity/refusal unit groups pass; all six Fourier,
nine forced-cycle and six harmonic groups pass, including the independent
small-force/large-unchanged-velocity refusal. The check measures consistency
at authored boundaries, not forward error, and does not infer unresolved
intra-segment excursions. Some nodal sampling configurations may still refuse.

## Complete spatial Maxwell histories on native rigid-body graphs

The spatial viscous-network owner now exposes network::relaxation.
Each RelaxationNetworkLink retains the original finite-body geometry and
full Kelvin matrix and adds explicitly supplied spatial Maxwell branches.
A branch stores its original six Gram rows L, positive relaxation time tau,
and all six elastic energy coordinates a. Its law is
a' = L s - a/tau, wrench = L^T a, with storage |a|^2/2 and
heat rate |a|^2/tau. Original factor rows define the history basis.
Equal-time branches, rank-deficient factors and mechanically unobservable
histories remain distinct; no projection, rank cutoff or pole merging removes
their storage or heat. The existing strictly admitted spatial stiffness factor
can initialize a branch, and an explicit elastic displacement can initialize
its memory. No relaxed stiffness, rest geometry or fitted coefficient is inferred.

Original native wrench columns B and complete body kinetic coordinates y
produce one simultaneous system:
y' = B(f - C B^T y - sum L^T a), s = B^T y.
Identity storage, skew body-memory coupling, original full Kelvin resistance
and diagonal Maxwell rates enter the existing exact quadratic pHS owner.
Preparation caches actual port-slip and passive-wrench integrals directly;
integrated memory forces are not recovered by subtracting endpoint histories.
Each connector's original Kelvin and Maxwell positive source rows independently
integrate heat along that same complete trajectory. Signed actuator work remains
a separate observation. All integrated net impulses enter the original native
simultaneous network once, retaining full inertia, attachment transport and
the original world linear and orbital-plus-spin angular momentum gates.

A conservatively PROVED null initial branch stress L^T a0 = 0 follows its
exact homogeneous decay independently of the driven graph response, which lies
in range(L). The original TwoProd/TwoSum exactness-witness rules, including
protection against an underflowed product remainder, admit this separation;
a rounded zero or numerical tolerance cannot. Its complete six amplitudes and
independent positive heat are restored using cached decay/expm1 coefficients.
The full original graph operator and every state remain prepared. Failed
proof leaves the original dense evolution unchanged. This exact superposition
prevents large unobservable histories from contaminating a tiny history through
independently rounded dense transition columns. Cross storage and loss vanish
by the original constitutive equations, including under nonzero actuation.

Accepted fixed-geometry receipts expose complete body and memory endpoints,
per-link positive heat and signed drive work, memory storage change, native
kinetic change and independent native impulse work. Both native work accounts
and the complete source-coordinate storage change must balance the independently
integrated heat and actuation. Initial and final native port rates must agree
with the complete source solution. The public original-to-final memory-energy
helper retains (after-before)*(after+before)/2 before rounding; it does not
subtract rounded absolute stored energies. Zero duration preserves the complete
body/history bank and exact zero impulses and heat, including large held efforts.

The moving entry composes two complete source halves around the existing loaded
native body drift. Memory amplitudes retain their original first-body material
basis through that drift, with no material-energy change. Actual endpoint
geometry and world-held efforts are re-expressed at each source. This is an
explicit corotated rate model; it does not claim a finite-strain equilibrium
spring. The whole-step account independently retains original-to-final kinetic
and memory storage changes, applied load work, gravitational storage, internal
actuation, material heat and native drift defect. Accepted state and history
can be transferred together without copying through into_state_and_memory.
The caller must not reinterpret existing amplitudes under changed factors
without a separately admitted thermodynamic history transfer.

Shared native geometry preparation is extracted from the existing driven source
without changing its arithmetic. Source scratch and exact caches are reusable;
warm fixed-geometry calls evaluate contractions and one native impulse event,
with bounded linear output allocations and no exponentials or factorizations.
Complete body, connector, branch, byte and work bounds precede preparation.
Both twice the full body-plus-memory coordinate count and that count plus all
physical input coordinates and one must fit the existing 256-coordinate owner.
The moving budget includes both preparations and retained candidate histories.
Borrowed bodies, sources and original histories remain unchanged on refusal or
cancellation. Moving geometry still requires independent temporal refinement.
Material validity, temperature-dependent factor changes, equilibrium springs,
contact and caloric receivers remain separate explicit owners.

Verification: All six focused fixed-source groups passed natively after the null-history
correction: analytic single-arm recoil/storage/heat; original unobservable,
equal-time and tiny histories with driven superposition; zero-arm original
viscous parity; full original-coordinate tensor force/torque/history/heat/work
ODE and fixed-clock continuation; atomic memory/source/resource refusals;
and compensated original memory-energy change. Three moving groups passed
against the unchanged moving-source bytes: independent corotated Newton-Euler
and memory time refinement, loaded common-world-rotation objectivity, and whole
body/history cancellation and bounds. Nine pre-existing driven graph regressions
also passed after the shared native geometry extraction. Four hereditary finite
thermal integration groups passed against the corrected fixed source. These
are bounded isolated builds of the actual native MBD/time, pHS, linear-algebra
and material/thermal numerical sources, with no numerical substitutes; they are
not a full workspace build. Pinned Rust formatting passed for the new modules.

### Complete affine-force relative-speed bound

`IsotropicLinearDrag::affine_force_relative_speed_bound_m_s` admits the same
held bath, gravity, inertia and zero-couple inputs as the exact affine-force
flow. It supplies a conservative bound on the entire COM slip path, including
interior extrema that endpoint-speed checks can miss. With `lambda=c/m`,
`S0=|v0-U|`, and `A=max(|F0/m+g|,|F1/m+g|)`, positive exponential convolution
gives the monotone scalar comparison
`exp(-lambda*t)*S0 + A*integral_0^t exp(-lambda*u)du`; its maximum is the larger
of the initial and final comparison values. The implementation uses the
original compensated momentum subtraction and decay moments, including the
continuous zero-drag limit `S0+duration*A`. This is a floating-point model
upper bound, not an actual peak measurement or interval certificate. It does
not make the mechanics core a material-regime or force-bandwidth authority.
The focused affine regressions include a positive interior maximum between
two zero-slip endpoints, gravity cancellation, zero/tiny resistance, and
large common translation.

## Native angular history after a material-temperature update

Both rolling and torsion event states expose `continue_at_temperature` for the
original sealed source lifecycle continuation. Exact predecessor comparison
includes attached spring memory and all accumulated event/loss accounts. The
native state's original body-pair geometry, inertia binding and contact frame
are retained. Source and returned native state must be published together;
ordinary subsequent native advances continue to enforce held geometry.

This held-coordinate continuation changes no body momentum, pose, mechanical
clock or thermal store. The original source report supplies signed parameter
work. Cancellation is checked before and after constructing the immutable
successor; foreign predecessor history refuses. Validation in this increment
was source review, pinned formatting and recursive parsing, not an executed
cross-crate native test or a finite-rotation/thermoelastic validation.

### Held viscous kinetics on the original Maxwell graph

HeldRelaxationKinetics supplies one nonnegative Kelvin multiplier per original
connector and one positive relaxation-time ratio per original branch.
PreparedRelaxationNetwork::new_with_kinetics and the corresponding moving
entries apply these coefficients to the original reference laws. Each effective
time and its reciprocal must remain finite and positive. The existing
PassiveBushingMatrix::scaled owner admits each Kelvin form while retaining its
Gram directions and coordinate scales; no matrix or memory factor is refitted.

The effective Kelvin matrix supplies geometry preparation, full resistance,
positive source rows and integrated passive-wrench observations consistently.
The effective branch time supplies pHS resistance, independently positive
Maxwell heat and the cached exact null-history decay/heat. Every elastic Gram
row, original branch order and complete memory coordinate remains unchanged.
Changing only these dissipative coefficients changes no stored elastic energy
and introduces no parameter work. Original constructors and moving entries
delegate to the same implementation with identity kinetics.

The public allocation-free validate method checks the complete bank and all
effective Kelvin forms, times and rates without preparing dense flow. It shares
the preparation's coefficient admission helpers, so a thermal owner can admit
its final refreshed endpoint coefficients without advancing body, history or
clock. Prepared-source getters expose the actual effective Kelvin form and
per-branch relaxation time without allocation. Extra cached times remain
inside the existing complete native storage bound.

Both source halves of a moving native interval use the same borrowed kinetic
bank while retaining the existing full body/history/heat/work transaction.
The composing material and thermal owners supply subsequent banks from their
original calibrated reference laws; rates must not be compounded from earlier
shifted results. This native layer supplies no temperature calibration,
temperature-dependent stiffness, entropy model or history reinitialization.
Coefficient variation within an interval still requires independent temporal
refinement; exact held-source heat is not a thermal trajectory error bound.

Verification: All five focused kinetics groups passed natively against the final source:
bit-identical identity kinetics; nonuniform Kelvin and branch-clock scaling
against the independent full native-coordinate momentum/history/heat/work ODE;
changing clocks with original tiny null histories, zero-time storage invariance
and no compounded shifting; matching coefficient admission and atomic retry;
and both moving source halves with full tensors, native loads, mixed world/
material-held actuation and cancellation. All nine pre-existing fixed/moving
hereditary groups also passed (14/14 total). These bounded isolated runs use
the actual native MBD/time, pHS and linear-algebra numerical sources without
numerical substitutes; they are not a full workspace build. Pinned Rust
formatting passed for the two modified modules and two new test files.

## Forced response about an affine-constrained equilibrium

`ConstrainedQuadraticModel::prepare_forced_response` reuses the original admitted equilibrium, allowed-coordinate map, native joint mass, and original allowed-direction stiffness and velocity tangents. It prepares the original general quadratic forced-response factor for peak phasors with the `exp(-i omega t)` convention, including DC. Physical joint forcing is reduced by `N^T`; physical displacement is lifted by `N`. Signed stiffness and nonsymmetric or gyroscopic velocity terms remain unchanged. This local particular response does not establish transient decay, nonlinear stability, or finite-amplitude joint-limit reachability.

Reaction increments are recovered using the original QR owner with the admitted physical row and coordinate scaling, then replayed through the original `A^T` action in original row order. Publication separately checks the complete physical force balance, original affine position and velocity tolerances, and the caller's explicit bound on average ideal-reaction power. Reported reactions are increments about the mean equilibrium, and the signed input, inertia, velocity, stiffness, and reaction powers refer to the actual published physical phasors. Constraint multipliers never receive artificial mass or damping.

A fully locked mechanism has a legitimate empty reduced response and recovers support forces without probing fictitious dynamic directions. Repeated loads reuse the native dynamic and reaction factors. Singular resonance, unrepresentable arithmetic, exceeded work or payload bounds, violated physical residuals, and cancellation refuse publication without changing the source model. The existing modal and unconstrained response APIs retain their admission rules.

Five original-mechanism regressions and one additional native modal power regression were authored and independently source-reviewed. They remain uncompiled and unexecuted during the executor outage.

## Native four-coordinate adhesive temperature continuation

`jkr_combined::temperature` retains the original normal/vector-shear/twist
native pair through an authenticated material-temperature update. The existing
native history check admits the actual body association, held poses, marker
datum and all four rates first. The original combined source then re-evaluates
one relaxed contact radius and supplies all four signed storage changes.
No normal-only or independently refreshed mode response substitutes for it.

Bodies, momenta, mass/inertia, full 4x4 mobility, contact frame, material and
reaction arms, coordinates, rates, pending COM/angular travel, time, steps and
budgets remain fixed. Only the source binding and cached combined response
change. Native energy tolerance applies to parameter-work closure, and final
cancellation leaves no successor. A temperature update cannot repair stale
geometry or consume pending finite travel.

Three integration groups exercise all four loads and independent temperatures,
continued native recoil, full off-diagonal mobility, common native motion and
nonzero pending travel/chronology, exact round-trip retention, and transactional
foreign/wrong-source/stale-body/stale-geometry/cancellation refusal. Pinned
formatting, recursive parsing and source review passed. These cross-crate tests
were not compiled or run under the session's memory/disk constraint.

## Exact simultaneous affine force and uniform bath velocity

`IsotropicLinearDrag::step_affine_force_and_flow` integrates
`m dv/dt = F(t) + m g - c(v-U(t))` with both force and uniform bath velocity
prescribed by their interval endpoints. Mass, isotropic resistance, gravity
and the free-spin owner remain fixed. The prepared bath supplies the starting
velocity; continuation must prepare the next owner at the accepted ending
velocity. This is the native moving-damper law, not a constitutive assertion
about accelerating-fluid pressure, added mass or Basset history.

The original compensated small-exposure and stiff exponential bases now also
integrate the relative forcing `F/m+g-dU/dt`. Positive-square slip moments
supply heat independently of storage. Applied work integrates `F(t) dot v(t)`;
signed bath work integrates `-c(v-U(t)) dot U(t)` with the actual velocity
history. Pairing drag impulse with the endpoint-mean bath is not substituted
for that source integral. Momentum, represented kinetic/gravity storage,
free-spin decay and the work residual keep their original owners.

The accompanying whole-path bound uses the endpoint norm envelope of
`F/m+g-dU/dt` in the original positive exponential convolution. It covers
interior slip peaks even when both endpoint slips vanish. It is a conservative
upper bound, not a reported peak or a physical bandwidth certificate. A held
bath delegates to the original affine-force entry; zero drag delegates before
forming a bath difference, preserving exact ballistic independence.

Focused regression fixtures compare both analytic bases with an independent
absolute-velocity ODE and power quadrature, compose actual force/bath
endpoints, check Galilean heat/impulse/work transformations, resolve tiny-drag
heat, retain exact forced co-acceleration, and exercise interior-slip admission.

### Complete lossless equilibrium on the spatial material-memory graph

SpatialMemoryEvolution distinguishes a finite positive Maxwell time from an
explicit Equilibrium channel. SpatialRelaxationBranch::from_equilibrium_gram_rows
retains all six authored factor rows L and complete elastic amplitudes a. The
same native work port gives a'=L*s, wrench=L^T*a and storage=|a|^2/2; its resistance
and positive heat rows are exactly zero. No infinite/large surrogate time,
rest displacement, Log(R) spring or material coefficient is inferred. Initial
equilibrium memory is supplied explicitly, and moving geometry retains the
existing objective corotated rate-law semantics.

Equilibrium and Maxwell channels share the original simultaneous full-graph
pHS flow, held drive, native wrench transport and independently gated kinetic,
memory, work and heat accounts. Every channel remains in the complete state
and resource bounds. Exactly witnessed null-stress equilibrium histories
superpose unchanged, including beside tiny relaxing histories; they generate
no heat and are never dropped or projected into a replacement basis.

The typed evolution getter retains the physical distinction. The original
relaxation_time_s getter now returns Some(finite_time) for Maxwell and None
for equilibrium. Prepared effective_evolution distinguishes a lossless channel
from an absent index; effective_relaxation_time_s returns only finite Maxwell
times. HeldRelaxationKinetics still requires a finite positive time ratio for
each original channel, but leaves Equilibrium exactly lossless. Kelvin scales
and Maxwell times retain their prior coefficient admission. A kinetic update
does not change equilibrium factors, stored energy or complete history.

PassiveBushingMatrix::from_gram_rows admits an explicitly supplied physical
factor without a second semidefinite factorization. The compensated full
L^T*L supplies its native SI matrix; original physical rows and coordinate
scales remain retained, including dependent/zero directions and weak positive
rows. Finite coefficients, physical congruence and representability are
checked in fixed-size allocation-free work. Existing matrix constructors keep
their original arithmetic. Scaling a factor-authored form preserves those
same physical directions rather than refactoring its rounded matrix.

Verification: Four focused additions passed natively: exact lossless forced
oscillation with signed work and zero heat; complete equilibrium/Maxwell/Kelvin
recoil against an independent forced ODE and exact-clock continuation; strong
lossless null history beside tiny Maxwell memory/heat with shifted kinetics
and atomic retry; and unchanged original Kelvin factor rows under physical
scaling, including dependent and weak directions. All fourteen pre-existing
fixed/moving hereditary and kinetics groups also passed (18/18 total). These
bounded isolated builds execute the actual native MBD/time, pHS and linear
algebra numerical owners with no numerical substitutes; they are not a full
workspace build. Existing thermal getter assertions were adapted to the typed
finite-time API without changing thermal mechanics or heat routing.

## Actual native motion with retained combined adhesive history

`PreparedRigidJkrCombinedPair::continue_axial_geometry` turns accepted COM
travel and commuting axial spin into actual poses. Its shared
`continue_at_geometry` core also admits caller-supplied compatible poses and
a tangent gauge change. Both use the original material markers and contact
datum, rebuild the actual complete mixed-unit mobility, retain the complete
source history with objective shear transport and unchanged chronology, and consume pending travel once.
The existing ordered geometry capture is shared with ordinary preparation;
noncontiguous original assembly indices survive pair-local recapture.

Eligibility is historical: each body's angular velocity and world spin
momentum must stay on the original normal throughout every accepted held
phase. Later cancellation of transverse excursions does not repair that path.
Endpoint attitudes must match the sealed integrated exponential-map rotation
within representation roundoff. Common material spin is the mean of the two
retained signed axial angles, independently of tangent gauge selection.

Double-double audits independently check actual COM displacement, linear and
world-spin momenta, total orbital-plus-spin angular momentum, kinetic energy
and normal approach. The original source projection continuation checks every
stored-energy component. Projected kinetic repartition remains diagnostic;
no impulse, gap correction, geometry work or heat is inferred.

Four new integration groups specify 48 successive free normal/twist phases,
a fully loaded shear/twist contact with actual common spin and changed tangent
gauge, original noncontiguous body indices, false/cancelled endpoint refusal
with retry, and cancelling transverse angular excursions. Pinned formatting,
recursive Rust parsing and scoped whitespace checks pass. These integration
groups remain uncompiled and unexecuted under the workspace resource limits.

The normal remains fixed and original transverse marker-collinearity admission
is retained. This owner does not provide arbitrary finite attitudes, rolling,
material-marker migration, detached trajectories or moving thermal faces.
