# CONTRACT: fs-mbd

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
rigid-rotor and transient owners described here. Finite internal shaft actuation
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
physical operators. Cancellation or failed balance publishes no partial state.

The transient tests exercise full-mass impulse-to-trajectory response, native
rigid-rotor angular momentum and energy, analytic whirl and damped ring-down
refinement, nonsymmetric active work, and singular/underresolved/cancelled refusal.

`rotordynamics::rigid::RigidRotor` supplies a concrete physical assembly. It
lowers the original full centroidal inertia by the supplied proper reference
rotation and admits axisymmetry under an explicit tight tolerance. Its four
stationary coordinates are lateral translation and tilt; finite bearings use
their actual axial stations through `A^T B A`, retaining force moments.
Signed-speed unbalance uses `A^T F`, and recovered bearing work keeps its sign.
This is a rigid, constant-speed linearization, without flexible-shaft,
spin-acceleration or asymmetric-rotor Floquet claims.

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

`gravity::CentralGravity` owns frame-checked Newtonian acceleration and specific
potential, with caller-supplied source position and G M. Rocket environments
reuse this owner. `UniformGravity` binds the existing acceleration to an explicit
potential datum. Neither infers Earth, a radius, softening or a rotating frame.

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

Pair work is admitted before allocation. Cancellation is polled per body/pair and
at publication; refused allocation, arithmetic, frame, budget or separation leaves
the immutable input system untouched and returns no candidate. Both endpoints and
the complete straight relative drift segments must clear the separation cutoff,
preventing a split drift from tunneling through a source. This is not a validated
enclosure of a true curved encounter: callers must resolve near-field motion.
`tests/mutual_gravity.rs` covers analytic finite-source recoil, independent force/
potential derivatives, the two-body Kepler period and second-order apocenter error,
20-orbit energy/momentum behavior, frame/boost/permutation covariance, reversal,
native rotated-inertia spin, hidden interior source crossing and late cancellation
(G1/G2/G3). This monopole lane does not include finite-body effects. Self-energy
within a body, collisions, deformation, external loads, variable mass, adaptivity
and FMM remain outside this closed-system map.

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

## Configuration-dependent loads on the original rigid integrator

`RigidBodyIntegrator::step_verlet_with` composes two endpoint force/torque
half-kicks around the original momentum-preserving full-inertia free map.
The first load belongs to the supplied state; the callback evaluates the
actual drifted endpoint and returns its original load receipt with its
wrench. Gravity acts once through the native owner. Each body torque uses
its own endpoint axes, so its world impulse is rotated at that endpoint.

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
The single-plane impact and shared normal/frictional/stateful contact mechanics
consume these increments rather than subtracting whole-tree totals. Native
velocity updates, endpoint speed limits, constitutive tolerances and existing
roundoff allowances remain unchanged. Focused regressions cover boosted contact
loss, offset-world angular impulse, internal joint work and recoil, and an
unrepresentable requested velocity increment.

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
by large bulk motion. Signed parameter work is evaluated in that same frame
with the actual linear-momentum-residual correction; it remains available when
rounded lab kinetic totals lose a small rotational-energy change. Positive
and negative work both consume the absolute budget. No heat partition follows.

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
momentum residuals independently audit the entire step. Lab-frame kinetic
subtraction can lose small signals under a large common boost; the separate
conditioned parameter-work receipts retain their own numerical meaning.

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

Caller effort is evaluated before reactions are added, at native elapsed times
zero and half-duration. The first callback sees the virtual constrained rates
under midpoint inertia. Independent stage powers keep caller actuator work
separate from reaction work; the retained original native actuator-work receipt
continues to include their sum. Original generalized-effort limits admit the
combined actuator and constraint load.

Separate complete-window caps count absolute parameter and reaction work.
The last transfer receives only each remaining allowance, including reaction
work already supplied during smooth motion. Opposite signs cannot hide work.
Both transfer receipts, both stage reactions and the native motion result are
retained, with independently calculated final energy and world-momentum residuals.
A failed stage/helper/native step or cancellation returns no state; typed helper
errors survive the callback adapter. Callback-owned side effects cannot roll back.

Second-order refinement requires smooth resolved loads and prescribed tensors,
an O(duration²)-accurate midpoint sample and adequate mode resolution. This
does not turn native explicit midpoint into an energy-exact or reversible method.
The focused G1 tests compare driven nonuniform inertia against scalar reduced
gear equations and analytic reaction torques, then state-dependent damping
against closed-form motion and separately integrated work. Both use 8/16/32
steps over one second; empty-row equivalence, gross work and midpoint cancellation
are also exercised. This remains a constant ideal-joint model, excluding
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
at every observed checkpoint. This supplies smooth regular bilateral distance
graphs; spherical/weld graphs, redundant constraints, unilateral cables,
contact/impact coupling, moving lengths, exact energy conservation,
symplecticity and validated whole-path enclosures remain separate work.

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

During support, only the normal applied force is balanced. With tangent inverse mass `kappa = 1/m + radius^2/I`, rolling requires traction `-(F_t/m)/kappa` inside the circular friction limit. Collinear sliding retains analytic constant-force motion, exact slip-arrest times and resolved reversal. Frictionless support admits arbitrary tangent drive.

Resolved changing-direction sliding uses the original isotropic Coulomb disk in a velocity-implicit discretization, bounded by `max_sliding_step_s` and `max_sliding_substeps`. Actual translation, world momentum, contact-point slip and torque include the original spin inertia. Each accepted step uses its resolved constant force in the native trajectory; endpoint Coulomb feasibility does not assert an exact continuous force direction. The curved trajectory is first order and arrest is resolved at a timestep boundary. Refinement is necessary. Coarse trials that would make friction assist the initial slip are halved, including apparent endpoint arrest. Supercritical transverse drive cannot be silently converted to persistent rolling. Rejected trials count toward the work cap.

For distinct coefficients, the transverse noncapture condition retains the kinetic limit: driving at or above that limit cannot acquire a finite-time zero-slip endpoint merely because static friction could hold a separately supplied resting state. Thus the hysteresis band admits both persistent rolling and continuing slip according to actual history. The exact critical transverse case also retains its physical noncapture condition. The five focused groups in `tests/sphere_plane/support_friction.rs` check analytic breakaway, arrest/re-stick/reversal and zero-kinetic capture; unchanged impact response; curved noncapture; and moving-support work, split continuation and atomic refusals.

All five and the fifteen prior sphere/plane groups pass with the pinned compiler against the original selected native dynamics/contact sources and full real tribology/math dependencies. This is focused native verification; it does not establish the full Cargo/DSR workspace gate.

Both routes evaluate lost work from actual midpoint relative slip on each constant-force segment, never from the fictitious impulse-law loss at the force-driven free velocity. Curved-step work belongs to the declared time discretization and converges with its sliding timestep. Rolling traction has zero relative loss; the moving plane supplies signed work. No mechanical loss is automatically assigned to heat.

Attitude uses an exact exponential for commuting spin/torque and bounded second-order Lie midpoint substeps otherwise. Newly supplied penetration refuses. Each curved supported segment checks its complete quadratic normal-gap path against the original floating geometry allowance, in addition to actual endpoint kinematics. Private continuation retains only bounded contact roundoff; no pose or momentum is projected.

Every accepted advance independently checks total linear impulse, orbital-plus-spin angular impulse about the world origin, and `K_after - K_before = external_force_work + plane_work - dissipated_energy`. Cancellation, unresolved geometry/work, nonadvancing clocks and exhausted event, attitude or sliding budgets publish no partial state. Accumulating rebounds do not alter restitution.

Fifteen focused native integration groups pass: the eleven original ballistic, support, collinear, incline, spin, frame, attitude and refusal groups plus four curved-flow groups. An independent separable continuum solution for `u' = a - B*u/|u|` verifies velocity, position, spin and loss convergence at 16/32/64 steps, with error ratios 1.8–2.2; a separate 32/64/128 ladder verifies arrest followed by driven rolling. Coarse-step passivity, supercritical non-arrest, moving-plane work, split continuation, cancellation and late sliding-budget rollback are also exercised. These are focused actual-source native tests, not a full Cargo/DSR workspace gate. General multicontact Moreau-Jean dynamics, external torque, anisotropic bodies, multiple planes, finite-mass plane recoil, rolling/torsional resistance, arbitrary-shape certified CCD, source calibration and contact adjoints remain outside this owner.

## Native rigid-body simultaneous frictional impacts

`frictional_impact::apply_rigid_frictional_impacts` binds the simultaneous
normal/Coulomb owner to original rigid-body states. Actual masses, full body
inertias, orientations and COM-relative contact arms determine the full
normal/tangent Delassus matrix. Assembly contracts only shared-body influence
blocks, avoiding products over unrelated body coordinates. Stable contact IDs
also canonicalize native impulse accumulation.

Dynamic pairs receive equal/opposite original point impulses; translating
prescribed planes supply an explicit external impulse, angular impulse and
signed support work. Independently reconstructed native endpoint velocities,
kinetic/midpoint work, total linear/angular momentum and whole-event energy
must agree before any body is published. Initial penetration, unresolved
contact or any late balance failure leaves all bodies unchanged. Distinct
dynamic witnesses are not silently moved to a fabricated common point; a
tangential couple due to an admitted finite gap can still fail angular balance.

Four native physical tests cover off-center rotated anisotropic inertia with
independently derived normal/tangent cross-coupling, sticking/restitution and
moving-support work, a simultaneous three-body chain with exact input-order
invariance, and penetration/redundancy/cancellation/late angular refusal.
The contact core's SPD and determined-response scope remains in force.
This point-impact phase holds poses fixed, requires supplied geometry and
coefficients, and makes no automatic finite-patch preload, torsional-friction,
collision-discovery or impact-to-heat allocation claim.

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


### Native redundant support and prescribed-support work

`frictional_impact::apply_rigid_frictional_impacts` consumes the contact core's explicit impulse/velocity determinacy receipt. Its original shared-body Gram mobility makes a unique contact velocity image sufficient for unique native body momenta, even when individual contact impulses are indeterminate. Such impulses remain an admissible representative; the report does not promote them to unique local pressure or heat.

A velocity-only event carries the contact core's explicit impulse-work range witness, proving that its total relative work is independent of the selected contact-impulse representative. The original Gram mobility then also determines total prescribed-support work from the actual kinetic change. This admits compatible translating or rotating supports, mixed finite-body contacts and independent islands without requiring a common support translation. Exactly frictionless tangent coordinates can be excluded from the work proof, so distinct harmless tangent belt speeds remain admissible. Incompatible redundant source work returns `Contact(UnresolvedWork)`; the directional owner retains the original transformed error. Individual support loads and source-work allocations remain representatives. Existing actual endpoint-velocity, kinetic/work, linear/angular-momentum and cancellation gates still precede atomic body publication.

Native G1/G3 groups in `tests/frictional_impact/redundant.rs` cover a finite-height, three-point body in stick and gross slide with independently derived load transfer, energy and momenta; frictionless duplicate normal rows with restitution; conveyor work; and a rotating platform coupled to a finite-body collision under frame changes and contact permutations. Directional regressions retain the transformed witness and circular-limit behavior. Generic rank/feedback, configuration, contact discovery and physical heat-partition boundaries remain explicit.


### Automatic native sphere collision trajectories

`frictional_impact::spheres::advance_sphere_contacts` advances complete COM-centred native spheres under explicit uniform `Gravity`. `advance_forced_sphere_contacts` additionally admits one finite constant world force at each original centre of mass, in body-slice order. Every pair has an immutable coefficient/restitution law; stable sphere IDs canonicalize discovery and simultaneous impulse order. Equal force/mass accelerations retain the original straight relative-path query. Different accelerations use `fs-contact::swept_sphere::accelerated::sphere_path` to isolate every interior distance extremum and bracket first entry on each complete quadratic relative path.

No pose or witness is projected. Up to 64 outward floating-time neighbours resolve nonnegative actual witness gaps AND complete path gaps; unresolved geometry refuses. Original radial material arms retain mass, principal inertia and orientation. Nearby contacts enter one simultaneous normal/Coulomb event at the declared spatial resolution. Native free rotation keeps its original approximation and acceptance gates.

Gravity and additional centre-force work use actual accepted displacements. Linear and world angular impulses remain separate. Gravity torque integrates the original total mass moment, including total applied centre force; each body-specific force torque integrates its actual piecewise quadratic path across collision-induced velocity changes. Each segment and the full trajectory close kinetic work, linear momentum and orbital-plus-spin angular momentum before atomic publication. Impact loss remains unallocated mechanical dissipation.

Four new G1/G3/G4 groups extend ten original collision/gravity cases: independently forced elastic entry with analytic gravity/force work and angular impulses; manufactured curved 3-D frictional entry with an independent original-inertia Coulomb impulse and post-impact trajectory; accelerated departure/return, permutation and split continuation; and late event-budget/cancellation rollback, deterministic retry and refusal of unresolved loaded resting pairs. These focused original-source native tests do not establish a Cargo/DSR workspace gate.

This bounded floating analytic lane retains all-pairs discovery. It does not claim arbitrary shapes, interval-certified CCD, force fields, time-varying forces, imposed spin torque, persistent supports, rolling resistance, thermal allocation or full Moreau-Jean dynamics. Zero-speed contact requiring a sustained reaction refuses rather than integrating into overlap. The broader `frankensim-ext-contact-nonsmooth-lane-oh0i` remains open.


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
bounded spatial broad phase with the original native force/contact/drift step;
it does not introduce a second impulse solver. Each body retains its actual
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

The caller's maximum step controls the original first-order force/contact/drift
approximation. Geometry retries do not establish an adaptive time-error
estimate. This numerical lane supports straight drift of spheres and fixed
normal translating planes; it does not claim interval-certified CCD, arbitrary
shapes, rotating supports, deformability, adhesion, rolling/torsional resistance
or thermal allocation. The original simultaneous solver's coupled-contact cap
still applies. The broader contact and multibody beads remain open.

Six focused native G1/G3/G4 groups cover resting support with nonzero material
restitution, falling contact and settling, two-layer stack load transfer,
independently derived incline rolling/sliding trajectories, conveyor impulse
and work, hidden interior collisions with encounter-time retries, and complete
rollback after retry/plane-query exhaustion or cancellation. All six pass with
the original native dynamics, contact, Newton/Krylov, geometry and math owners;
the seven shared broad-phase groups pass separately. This selected-source
native verification uses the pinned nightly. It does not establish a full
Cargo/DSR workspace gate; the required sibling constellation was unavailable.

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
gates. The native static capacity bounds the required traction at initial,
midpoint and endpoint states. An infeasible reaction refuses instead of forcing
rolling or substituting an independently saturated friction component.

The material coefficient remains owned by the original card. Each analytic
two-axis midpoint solve uses its native moment and native heat. The user sets
the time-resolution ceiling; an additional derived limit
`c * lambda_max * h <= 1` keeps every rolling eigenmode's decay factor
nonnegative. Work is bounded by an explicit substep cap. Viscous arrest remains
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

This is one held-geometry, held-normal-load pure-rolling phase. The surrounding
owner still advances poses and normal contact. It does not solve contact
discovery, lift-off, externally driven no-slip constraints, simultaneous
microslip/spin, multiple patches, or arbitrary-speed viscoelastic rolling. A
second shear reaction cannot be superposed on this same phase. Frozen-geometry
time accuracy is second order; this does not certify a complete trajectory.
G1/G3/G4 groups in `tests/rolling_friction.rs` use independently eliminated
principal-direction mechanics for both-body recoil, a three-rung exponential
motion/heat refinement ladder, nonreversing decay, exact zero limits, world and
surface-order covariance, source/traction/mode refusal, and loaded checkpoint,
late audit, cancellation and exact retry.
A seventh group requires exact isolated-pair recoil, heat and history despite
an energetic untouched body and nondefault contact indices.
All seven groups pass in the focused native original-source build. The independent
8/16/32-step ladders require error-reduction factors in `[3.95, 4.05]` for both
motion and heat. This is separate from a full Cargo/DSR workspace gate.

## Material-derived supported sphere rolling

`sphere_plane::rolling` advances the original supported sphere pose and momenta
under the existing material-derived Kelvin–Voigt rolling moment. The original
Hertz patch, source density/radius, held temperature and frequency, and
quasistatic validity remain binding. The constant-load vector equation
`(m + I/R²) du/dt = F_t - (A N/R) u` is integrated analytically, including
driven turning and reversal. Both endpoint static-traction norms bound the
entire exponential path; infeasible rolling refuses rather than projecting
away slip. Native linear and world-origin angular impulses, actual force work,
moving-plane work and continuous bulk torque work are audited. The nonnegative
heat integral uses exponential moments, not an energy defect. Normal spin is
preserved. Changing-axis attitude alone uses bounded Lie midpoint steps;
commuting attitude is exact. This is a homogeneous-sphere, held-load no-slip
phase with asymptotic coast-down, not sliding, impact, finite-time viscous
arrest, torsional loss or thermal transport. Generated heat belongs to the
card's deforming sphere and is already included in the native advance loss.
Focused physical groups cover exact coast-down, driven work, world/boost
covariance, continuation, elastic limits, source/traction refusal, late
cancellation and attitude refinement.

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
The rocket staging adapter keeps its dry/propellant partition and inertia
assembly admission and uses this same release kernel with stack workspaces.
It reconstructs receipts after conversion to its stored velocities; the
two-piece rocket path remains allocation-free.

Seven native assembly cases check an independent full parallel-axis assembly,
zero-impulse rigid motion inheritance, shared-body impulse work against native
sequential impulses, closed-form orbital-to-spin latching, frame covariance,
small loss under large translation, and atomic limits/cancellation/retry.
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
