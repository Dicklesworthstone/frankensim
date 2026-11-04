# fs-flyer — CONTRACT

Bead: frankensim-wf-root-guzez.4.1 (E3.1, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md §5.1
(ROUND 6 steady state). Evidence: data/wright-flyer/{flyer-reference,
canard-mechanics,geometry-conventions,frame-conventions}.

## Purpose and layer

L4. The aircraft-level Wright Flyer model. E3.1 ships the `FlyerDesign`
schema (Round-2 canard/mechanics/structure fields), control-topology enums,
admission with typed refusals, the component mass/inertia build-up
(cross-checked against the published single-lineage inertias), and the
derived-quantity panel. The later E5 lifecycle engine now also lives here;
the browser boundary remains in `fs-flyer-wasm`.

## Public types and semantics

`fieldsvc` evaluates its bound horseshoe and flat-ground images through the
same Rosenhead–Moore finite-segment kernel as `fs-vpm`. Velocity gradients,
vorticity, strain, Q and lambda2 use exact spatial derivatives; the independent
finite-difference divergence overlay remains a diagnostic with its declared
step. Geometry and core scaled by s give velocity /s and gradient/vorticity
/s² at fixed circulation. Positive-core axis and endpoint derivatives remain
finite and retain the core's nonzero vorticity. Unregularized on-segment
derivatives are masked invalid. The core reliability mask covers every
selected span, trailing leg and image segment; it does not claim compact
support of vorticity. Invalid core lengths and non-finite bound inputs refuse.
The field digest domain is v2 for this intentional physical kernel change.
`tests/fieldsvc_battery.rs` checks scaling, independent velocity differences,
axis/endpoint behavior, image cancellation and selected-source masks (G0/G1).
`fieldsvc::PreparedWakeField` additionally binds the actual connected near/far
`fs-vpm` sources to this same return set. Preparation stores physical segments
once, using each source's core and elapsed diffusion; subsequent queries reuse
that buffer. Images across the field service's z=0 plane use polar velocity and
axial vorticity reflection without copying sources or entering physical ledgers.
All retained physical sources must stay above the plane when images are active.
The connected sheet includes its bound closing line, so the legacy payload's
prescribed horseshoe cannot be selected alongside it. The source snapshot hash
binds both the original field state and the prepared near/far/image identity.

Sampling uses exact retained-source derivatives and independent velocity-only
finite differences of that same field. Its cost is O(points * retained edges),
with no multipole-derivative error claim or recovery of pruned edges. Explicit
source-count and source-visit caps precede allocation; the sample's eight-visit
bound includes center derivatives, six velocity probes and core masking for
each selected physical/image system. Queries poll cancellation every 256 source
visits and at point boundaries, and publish no partial field. Core masks use
each segment's actual radius. Active wake images are required for a total-flow
claim, along with the physical wake itself.

`tests/fieldsvc_wake.rs` exercises actual near/far partition invariance, full
Jacobian velocity differences, plane/reflection identities, bound-line
ownership, finite-core axis curl, singular masks, budgets, cancellation/replay
and immutable source snapshots. The old ring payload still has no physical
wake history; callers must supply it explicitly. Neither source path claims a
propeller field, volume pressure, resolved vortex-tube dynamics, or completeness
beyond its declared source components.

- `FlyerDesign { wing, canard, rudder, lateral, warp, masses… }` — frd-body-v1
  coordinates, +x forward from the wing leading edge, SI units, radians.
- `LateralControlTopology::{WarpWithSlavedRudder{ratio}, WarpIndependentRudder}`;
  `WarpStructureMode::{FlexibleTrussed, Rigidized}`.
- `MassComponent { mass, position, gyration_sq, body_from_component }` —
  six-bucket engine convention included; squared principal gyration radii
  and an explicit `So3` frame carry each component's intrinsic tensor.
- `mass_build_up()` → gross mass, gross CG and full body inertia about the
  gross CG. `cross_inertia_kgm2` retains rotated intrinsic products and the
  negative parallel-axis products;
  `inertia_tensor()` reconstructs the symmetric matrix and `rigid_body()`
  prepares the original six-DOF body directly from it.
- `derived_panel(hinge_ratio)` → canard volume, naive two-surface neutral
  point, fixed AND free-control static margins, thin-plate hinge-moment
  gradient 2π·(x_h − ¼).
- `digest()` — canonical content identity under
  `org.frankensim.fs-flyer.design.v1` for aligned components; explicitly
  oriented components select `design.v2` and bind all canonical component
  quaternions (a PhysicalScenarioId ingredient).
- `FlyerDesign::reference_1903()` — the dossier-valued reference config.
- `simloop::SimLoop` / `SimStateOut` — 120 Hz lifecycle state including
  longitudinal state plus Estimated reduced roll/yaw state. Snapshot v2 is a
  14-float layout whose first 12 words are the unchanged v1 prefix and whose
  appended words are roll attitude and heading.
- `spine::step` / `spine::advance` — deterministic six-DoF fixed-step
  integration. The array boundary is admitted into canonical `fs-ga::So3`
  before either load callback runs; body angular velocity advances through the
  `fs-time` body/right Lie step, and returned quaternions use the canonical
  representative of the rotation double cover. Zero-horizon advance and replay
  still perform body/state/timestep/group admission and return that canonical
  state without invoking the load callback.
  The second load evaluation receives a first-order endpoint velocity predictor
  including both halves of the initial moment impulse; smooth angular drag and
  forces coupled to angular velocity therefore retain second-order convergence.
  The accepted rotational update delivers the trapezoidal world angular impulse
  `dt/2 (R_start tau_start + R_end tau_end)` to roundoff. This is a finite-step
  balance, not an exact solution for general varying loads. The spine battery
  checks exponential angular decay, its integrated angle and coupled translation,
  and a constant world-frame torque on an asymmetric rotating body (G0/G1).

`spine::RigidBody` owns the shared, immutable `InertiaTensor3`. `new(mass,
principal_moments)` explicitly declares principal body axes;
`from_inertia_tensor(mass, tensor)` retains arbitrary products of inertia.
Principal-frame preparation occurs once, outside all tick loops. Ordinary
flight torque kicks use the full inverse tensor, and native free rotation
transports angular momentum through the same fixed body-frame transformation.
Replay and partitioned flight carry that same body; no geometry or load axes
are relabelled. The diagonal path retains its original arithmetic and digest.
G1/G3 fixtures independently check component raw moments, the actual
design-to-flight path, body-coordinate covariance and cross-axis motor recoil.
`MassComponent::body_from_component` maps its declared principal axes into
the unchanged aircraft body coordinates. The existing `So3` and `Mat3` owners
rotate the component's diagonal intrinsic form before mass-weighting and
parallel-axis assembly. Positions already live in body axes and are not rotated
again. This cold assembly introduces no eigensolve or new dynamics owner.
Aligned components retain their original arithmetic, reference numbers and
v1 digest. Nonidentity frames use the v2 identity with canonical quaternions,
so opposite quaternion representatives cannot create different scenarios.
Idealized point/line components with zero principal moments remain admitted;
the gross body still goes through the original positive-tensor preparation.
Nonfinite derived gross moments refuse before publication.

G1/G3 design tests compare all nine tensor entries with independent volume
quadrature of an oblique homogeneous cuboid, retain the pilot's original mass
and inertia, and feed the result through an actual flight torque step. G0
checks cover canonical frame identity, degenerate component moments and
unrepresentable mass properties. The input principal moments remain declared
data or the documented reference calibration; a rotation does not validate
that calibration or supply missing geometry, moving masses or deformation.

## Section pressure drag through stall

The parasite ledger's spatial `wrench_at` evaluates the material-point velocity
relative to local air and returns force, its moment about the supplied origin,
and dissipated relative power. Co-moving flow yields zero force without a
direction singularity. Longitudinal loads sample incident air and pitch-point
velocity at the declared parasite line, retaining the accepted wrench for
distributed roll/yaw loads. The shared six-DOF propulsion owner samples world
wind there and includes all rigid angular rates. Its success path
does not allocate per-component reports or recompute a digest. G0/G3 tests check
spatial work, reversal, co-moving air and isolated local gusts through both real
aircraft load consumers. The constant scalar component coefficients remain an
Estimated isotropic model; directional areas, local propeller wash and resolved
interference are not inferred by this kinematic correction.

Longitudinal flight, distributed roll/yaw moments, and the native six-DOF
propulsion assembly consume the same effective strip incidence returned by
`fs-wing`. Steady coupling, evolving axial wakes and first-harmonic wakes carry
that incidence with their accepted circulation; it includes geometric twist and
trailing-sheet induction exactly once. Profile drag uses the existing section
law's smooth transition over `|alpha_eff|=0.30..0.45 rad`: the declared attached
`cd0_wing` blends to `Cd90*sin²(alpha_eff)`, the wind-axis drag projection of
`fs-airfoil`'s separated flat-plate normal force. Lift uses that same transition.
Fully separated flow therefore carries pressure drag instead of retaining the
attached drag floor. Induced drag remains the separately owned classical term.

The force opposes each strip's accepted chord-plane relative velocity; its
moment uses the same physical strip arm. Local relative work is nonpositive in
both incidence directions. Attached-regime arithmetic is unchanged. G0 tests
check plate force decomposition, even positive drag and continuous transitions,
then reconstruct drag and work in an actual separated aircraft load solve.
The attached coefficient, transition angles, constant quarter-chord moment and
flat-plate continuation remain Estimated model choices. This change supplies
the missing drag channel of that model; it does not establish historical
deep-stall validity, Reynolds-dependent viscous section data or resolved wake
pressure drag.

## Finite-strip rotational fluid inertia

`addedmass::assemble_analytic_strip` integrates each finite uniform strip's
normal translation, pitch about its explicit span axis, and normal-velocity
variation along its span. With chord `c`, span `b` and normal added mass
`m_a = rho*pi*c²*b/4`, its centroidal rotary terms are `m_a*c²/32` about the
span axis and `m_a*b²/12` about `span_axis cross normal`. The pitch coefficient
is the flat-plate limit of the ellipse coefficient in
[MIT 13.021, supplemental problem A1, page 3](https://web.mit.edu/13.021/13021_2004/suppnotes/supp04.pdf),
with ellipse semiaxis `a=c/2`; the span coefficient is the exact uniform second
moment. These are coefficients of the existing independent-strip potential-flow
approximation, not a full three-dimensional rectangular-plate solution.

`Strip` supplies orthonormal `normal`/`span_axis` directions plus the existing
centroid `control_gain` and explicit `control_angular_gain`. All three positive
kinetic-energy modes contribute to rigid, hinge and cross blocks together;
their shared added momentum supplies the original workless gyroscopic bias.
The joint solve also retains the rigid body's complete inertia tensor.
The actual extraction path uses these terms and records algorithm version
`config-resolved-strip-v2`; its existing wire schema is unchanged.

G0/G2/G3 tests in `addedmass_battery` cover the thin-ellipse pitch coefficient,
hinge/rigid recoil, finite-span subdivision, reference shifts, oblique frame
rotation, length-to-the-fifth rotary scaling, PSD, work closure and invalid
orientation. The intentionally changed joint and extraction goldens distinguish
the added physical storage. This entry point freezes geometry/control
Jacobians. Cross-strip fluid interaction, free-surface effects and
finite-aspect-ratio corrections remain outside the strip approximation.

### Finite moving hinges and configuration inertia

`addedmass::hinged::HingedStrip` prepares a physical body-fixed hinge line and
zero-angle strip geometry. Finite rotations move its centroid, normal and span
together; normal and angular control gains derive from that same geometry.
The reference strip must have inactive control fields, so no authored gain is
silently replaced. Multiple strips can share one active coordinate. The original
three positive kinetic modes and deterministic mass accumulator are reused.

`assemble_hinged_strips` returns the complete added bias for
`T(q,u) = u^T M(q) u/2`: body momentum transport, `-Mdot*u` on all rows, and
`dT/dq` on each control row. Analytic derivatives act directly on rank modes;
no derivative matrices, geometry vectors or runtime finite differences are
constructed. Bias power is `-u^T Mdot u/2`, which balances the changing kinetic
metric. A released hinge can have nonzero kinetic-gradient torque even when
its instantaneous rate is zero. The original joint solve consumes these same
body/control/cross blocks and bias without lagging their accelerations.

Each hinge is fixed in carrier coordinates; serial moving hinge chains, elastic
strip deformation and density time derivatives are not included. The G0/G1/G3
`addedmass_hinged` tests check finite-angle geometry and exact centrifugal
reaction, released-control torque, every bias component against independent
test-only mass differentiation, joint-solve energy rate, finite-strip
subdivision and original resource/geometry refusals. Fixed-strip goldens retain
their previous arithmetic and values.

### Fluid inertia on the native articulated tree

`addedmass::articulated::LinkedStrip` binds the original finite-strip geometry
to a native `fs-mbd::articulated` link. `assemble_articulated_strips` consumes
that owner's `GeneralizedKinematics`, shared with the complete dry mass and
bias extraction. Fixed/root links, serial revolute hinges, prismatic and
helical joints inherit the same native poses, local velocity Jacobians and
zero-generalized-acceleration bias. No second geometry or dynamics tree is
constructed. Authored strip control gains refuse rather than compete with the
native joints.

The original three positive local kinetic modes supply `J^T G J` and
`-J^T(G c + ad*_V G V)`, with `V = J u` and `c = Jdot u`. The Lie bracket
supplies both rotational and translational momentum transport. This includes
off-diagonal control inertia and the full analytic configuration bias, without
runtime finite differences or per-strip derivative matrices. The adapter alone
permutes native angular/linear columns into the existing Flyer linear/angular
block order. Applied loads remain zero in this fluid-only contribution; callers
combine the complete native dry metric and bias at the same configuration.

The original control/strip caps remain active. Density is constant in the
inertial bias; strip deformation and cross-surface fluid interaction remain
excluded. `addedmass_articulated` checks recovery of fixed/offset-link and
finite single-hinge results, two serial oblique hinges and mixed sliding/helical
motion against independent Lagrange bias differentiation, local strip energy,
nonzero control cross blocks, released-control torque, finite-strip subdivision
and resource/geometry refusals. These are idealized physics checks, not a
validation of historical geometry or aerodynamic interference.

### Simultaneous body and hinge motion

`partitioned::coupled::step` advances the existing six-DOF body and one active
hinge as a single seven-coordinate midpoint solve. Its callback receives the
midpoint pose, hinge state, body velocity and exact half-tick memory state, so
`assemble_hinged_strips` supplies the actual finite-angle mass and analytic bias
at each iteration. The original `solve_joint` retains all body/control cross
blocks. Applied loads have one explicit owner; duplicate nonacceleration loads
in the fluid assembly refuse, and the step supplies dry-body momentum transport.

Linear hinge stiffness and damping enter the same effective matrix as
`I_dry + c*dt/2 + k*dt^2/4`, avoiding an iteration restriction from those linear
terms alone. Nonlinear fluid, gyroscopic and load coupling must converge within
the declared velocity tolerances and bounded iteration count. Only an accepted
state receives the full exact memory transition. Refusal returns no partial
state. Orientation reconstruction uses the existing SO(3) exponential.

The dry model is a constant body tensor plus independent scalar hinge kinetic
storage. A moving dry linkage that changes the body/control mass matrix requires
the articulated owner; fluid coupling does not supply that missing dry geometry.
General nonlinear energy and world-momentum balance are refinement properties,
not exact finite-step guarantees. `partitioned_coupled` checks an independently
solved stiff-hinge recoil problem and damping work, second-order inertial-frame
translation and finite-hinge kinetic-energy convergence, exact memory timing,
and bounded atomic refusals.

### Uniform ambient flow in articulated inertia

`addedmass::articulated::assemble_articulated_strips_in_uniform_flow` extends
the original linked-strip assembly with a prescribed spatially uniform world
velocity `U` and its world time derivative `Udot`. In each link frame it uses
`W = V - (0, R^T U)` and
`c_relative = c + (0, omega cross R^T U - R^T Udot)`.
The complete bias is `-J^T (G c_relative + ad*_W G W)`; the existing rank modes,
native motion map and dry/fluid mass factorization retain their owners.
Uniform flow does not change the mass blocks. The still-fluid entry point and
the zero-flow branch preserve the original assembly arithmetic.

`UniformFlowAssembly` privately binds its load blocks, disturbance kinetic
energy `sum(W^T G W)/2`, affine generalized fluid momentum `sum(J^T G W)`
and work data to the same admitted native state. Its `reaction` method takes
the actual generalized acceleration in Flyer order and returns
`Q_fluid = bias - M_fluid * acceleration`, total world force, body power,
ambient source power `U dot F_world`, and an independently reconstructed
storage rate. Their balance is admitted against an arithmetic-resolution
tolerance; it is not a certified interval bound. Native root translation
columns recover world force correctly even with a nonidentity root mounting
transform. No reconstructed generalized velocity, extra solve, matrix
differentiation or second kinematics map is required.

The existing smooth load callback consumes the original returned load blocks
at its actual midpoint. Fluid density and strip geometry remain fixed.
Spatial shear, ambient fluid rotation, cross-surface hydrodynamic interaction
and density-change terms are outside this model. Pressure-gradient force on
displaced fluid remains a separately supplied external wrench; a displacement
volume is not inferred from the strip model. These disturbance-storage and
source-work diagnostics apply to smooth evolution and do not extend the
homogeneous energy receipts of the existing impact owners to moving fluid.

G0/G1/G3 tests cover zero-flow equivalence and comoving cancellation, the
independent pitching-plate normal force in a tangential current, rotating serial
linkages under Galilean and spatial-frame transformations, a nonidentity root
mount, independent material-point storage differentiation, and the actual
midpoint response to accelerating current with affine momentum and ambient
work closure. Reference formulations: Fossen's
[relative marine-craft equations](https://www.fossen.biz/html/marineCraftModel.html)
and MIT 13.012
[fluid-force notes](https://web.mit.edu/13.012/www/handouts/2003/Reading7.pdf);
the latter separate displaced-fluid pressure force from added mass.

### Moving dry links and fluid in one coupled step

`partitioned::articulated::step` advances the native free-base articulated
state, up to eight active scalar joints and 256 links, with exact linear-memory
transitions. Every midpoint prepares one `fs-mbd` generalized inertia/motion
map. The callback uses that same map for attached strips and applied link-origin
wrenches. `addedmass::solve_articulated_joint` sums the complete dry and fluid
matrices, retains every body/joint and joint/joint block, and applies each
owner's inertial bias once. Moving dry COMs and rotated intrinsic link tensors
therefore participate in the same recoil solve as the surrounding fluid.

The adapter explicitly permutes native `[angular, linear, joints]` coordinates
to the existing Flyer ordering. It reuses the original fixed-order Cholesky,
with a diagonal congruence for mixed units and refusal of unresolved normalized
pivots. It adds no fictitious diagonal inertia or numerical regularization.
The original constant-body solve retains its arithmetic and fluid goldens.

World-axis wrenches are moments about each actual link origin; native virtual
work projects them and explicit uniform gravity into generalized loads. Only
total internal effort, including declared passive impedance, is compared with
native joint effort limits. External forces do not consume that actuator bound.
Spring/damper terms are eliminated implicitly in the same mass solve; position,
rate and effort limits remain active. Each joint declares its own rate tolerance
in radians/s or metres/s. Bounded iteration and cancellation publish no partial
state. The pose update uses the original SE(3) exponential.

The `partitioned_articulated` tests compare the full adapter with native ABA
under gravity and link loads, two stiff dry/fluid rotors with an independent
Schur/Cayley solution and damping-work balance, and an oblique moving-COM chain
with independently computed strip energy/spatial momentum. Energy, linear
momentum and angular momentum errors all converge at second order in that free
trajectory. No exact finite-step nonlinear conservation, impact, joint-stop or
closed-loop constraint claim is made by this smooth coupled lane.

### Simultaneous contacts with complete articulated fluid inertia

`partitioned::articulated::impact::resolve` resets a caller-located event on one
native free articulated tree with up to eight joints, 256 links and 64 normal
contacts. Link-point/self-contact and prescribed moving-plane geometry come
from `fs-mbd::articulated::contact_geometry`, also used by the original dry
contact adapter. The same admitted link motion map projects every impulse row
and the complete linked-strip inertia. No scalar effective mass replaces the
coupled dry/fluid body and joint blocks.

`addedmass::PreparedArticulatedMass` owns one immutable, diagonally equilibrated
snapshot of the dry plus fluid matrix. Its original Cholesky factor is reused
for every column of `J M_total^-1 J^T` and the final accumulated generalized
impulse. Momentum uses the assembled matrix; positive kinetic energy uses its
factor. Configuration, fluid density or attached geometry changes require a new
snapshot. Force/bias fields do not act during an instantaneous reset, and
duplicated fluid applied-load fields refuse at this boundary.

The original `fs-contact` simultaneous normal law chooses all impulses together
and retains its global restitution/passivity admission. Equal/opposite internal
point impulses conserve combined spatial momentum. Prescribed planes supply
explicit external impulse and signed work. Actual returned native point speeds,
total dry/fluid linear and world-origin angular momentum, and complete kinetic
work are checked before publication. The energy change uses quadratic
polarization to avoid subtracting two large kinetic energies. Native joint-rate
limits still apply; an impulse is not an actuator force and does not consume the
continuous-effort limit. Cancellation, invalid geometry, penetrating witnesses,
singular active rows and unresolved endpoint balances publish no candidate.
Pose, joint positions and aerodynamic memory remain unchanged.

`articulated_impact` verifies simultaneous self/plane coupling with independently
computed dry and strip-local fluid momentum/energy, signed moving-support work,
the zero-fluid original ABA response, and a rebound law that must refuse global
energy creation. Atomic ownership/cancellation, exact gap admission and row
caps are also covered. This single-tree normal entry point does not locate collision times or
provide friction, persistent contact or penetration repair. The shared mechanics
also serves the multi-tree and Coulomb consumers below.

### Impacts between independently fluid-bearing articulated trees

`partitioned::articulated::impact::resolve_many` composes up to 128 native free
articulated trees, each within the existing eight-joint/256-link bound, at one
caller-located normal-impact event. Each tree's complete dry and attached-fluid
mass is prepared once from the same native map used for contact projection.
The assembly is block diagonal between trees before contact: cross-tree
hydrodynamic added-mass blocks are explicitly absent. Shared contact rows
retain cross-contact mobility through every common incident tree; no global
dense body-coordinate matrix or independent scalar effective masses replace
the coupled dynamics.

The shared point-impact mechanics serves both normal and Coulomb consumers.
It applies the accumulated impulse once per tree, checks actual native endpoint
rates, complete global kinetic work, prescribed-support work and world-origin
spatial momentum, and publishes all states atomically. Configuration and
aerodynamic memory remain unchanged. The existing single-tree `resolve`
delegates to this same owner; generalized travel-stop reset remains its existing
physical-mass bridge.

`articulated_impacts` compares two separate fluid-bearing articulated rotors
against independent strip inertia, joint response and orbital-momentum
formulas; a simultaneous three-body chain verifies shared-tree negative
cross-contact mobility and complete inelastic loss. System/contact permutation,
an uninvolved moving tree, late native-rate and fluid refusals, cancellation
and system caps are covered. The normal law retains its existing
singular-active-row refusal; contact detection, persistent contact and
hydrodynamic interaction between trees are separate capabilities.

### Coupled Coulomb impacts with complete articulated fluid inertia

`partitioned::articulated::impact::friction::resolve` applies one simultaneous
normal and two-tangent point-impact law across independent native articulated
trees. It reuses the original `ArticulatedFrictionalContact` and shared native
frame preparation, canonicalizes stable contact identities, and projects every
row through each tree's prepared complete dry-plus-fluid mass. One original
`fs-contact::frictional_multicontact::solve_frictional_contacts` call retains
normal/tangent, cross-contact, base/joint and dynamic-counterpart coupling.
Each tree receives one accumulated generalized impulse; no sequential tangent
clipping, mass regularization or continuous force/bias contribution enters
the reset.

The original law's unique-impulse versus unique-velocity determinacy receipt
remains public. Redundant reactions can establish physical motion without
establishing individual pressure or heat. For a velocity-only witness, active
prescribed supports must share one translation; independently moving supports
refuse without a work-determinacy witness. The common mechanics owner
reconstructs actual endpoint link velocities and audits complete dry/fluid
kinetic work, support work and spatial momentum before publishing all returned
states together. Actual rate limits remain enforced, while an impulse does not
consume a continuous actuator-effort limit. Pose, joint configuration,
aerodynamic memory and caller input remain unchanged.

`tests/articulated_friction.rs` supplies independent analytic oracles for an
off-centre fluid-loaded hinge with coupled normal/tangent recoil and
moving-support/frame covariance; equal-and-opposite friction between separate
fluid-loaded bodies with the zero-fluid native limit; a redundant four-foot
contact with determined motion and canonical identity order; and ambiguous
support-work, duplicate-identity and late-cancellation refusal. This bounded
event admits at most 21 contacts and an explicit scalar Coulomb coefficient.
It does not discover collision times, infer cross-tree hydrodynamic mass,
supply persistent contact or finite-patch spin, choose material coefficients,
or allocate impact loss to heat/deformation.

### Rigid linkages during complete-fluid Coulomb impacts

`partitioned::articulated::impact::friction::constrained::resolve` preserves
admitted distance, spherical point and weld constraints during one caller-located
event. Original native bilateral geometry and normalized Schur condensation
consume the augmented mobility from the same per-tree prepared dry-plus-fluid
mass. The original Coulomb law solves on the condensed mobility; all contact
and recovered bilateral impulses are accumulated before one update per tree.
No penalty stiffness, sequential velocity projection or dry-only reaction mass
is substituted.

Up to 21 contacts and 32 independent scalar bilateral rows (95 total) are
admitted, within the existing per-tree bounds, 128 trees and 4096 total links.
Uncancelled endpoint response scales use the complete mass and reuse each
endpoint solve when forming combined rows. Native actual rates, complete
kinetic/support work and spatial momentum still gate publication. Admitted
numerical bilateral work is reported separately and remains in the energy
residual; it is not heat. Reactions have N s units for point/range rows and
N m s for weld angular rows. The inherited native scalar velocity threshold
uses m/s or rad/s according to row kind. Empty constraints delegate to the
existing unconstrained consumer. The original unique-velocity verdict gives
no new authority to local contact or linkage loads.

`articulated_constrained` checks analytic fluid-loaded distance transmission,
distinct off-centre spherical/weld recoil, the objective geometric couple at
admitted closure error, simultaneous redundant contact/chain ordering and atomic
refusals. This adds instantaneous constrained impacts; continuous loop
integration, collision discovery and cross-tree hydrodynamic mass remain
separate.

### Sustained contact through the coupled fluid-aware midpoint

`partitioned::articulated::contact::step_with_contacts` advances one articulated
tree through sustained fixed-feature Coulomb contact using the original coupled
Lie midpoint integrator. Contact reactions use the complete dry-plus-fluid
effective inertia, including the original analytic fluid bias and implicit joint
spring/damper elimination. One simultaneous contact solve acts on half of the
effective contact mobility, because a full endpoint impulse changes midpoint
velocity by one half of its endpoint response. The reaction receipt is aligned
with the previous response that generated the accepted midpoint iterate;
actual midpoint normal, cone and complementarity conditions are reconstructed
independently before publication.

Resolved incoming closing active contacts require the separate physical impact
reset. Arithmetic-sized negative normal speeds remain unprojected and are
returned with per-contact allowances derived from uncancelled native point and
support motion and the original free velocity response. Solver tolerances and
absolute speed floors do not enter that admission. Unresolved friction
stopping/reversal, release/recontact, negative actual gaps and artificial
internal spatial couples refuse. Contact within a declared positive activation
band remains explicitly a finite-resolution reaction, not a repaired gap.

The required `ContactPath` callback must admit the complete original
constant-midpoint-twist / linear-joint configuration path, continued witness and
frozen-frame validity, and omitted-feature coverage. Initial, actual midpoint
and endpoint nonpenetration checks alone do not establish continuous collision
detection. The callback's admission and a final cancellation check precede
publication; inputs remain unchanged on refusal.

The result exposes actual signed reaction, normal, tangential and translating
support work, plus the external spatial impulse about the world origin. The
auxiliary Coulomb solver's impact-loss diagnostic is not physical friction heat.
Redundant reactions retain the original unique-velocity versus unique-impulse
distinction and the common-support-translation restriction for determined work.
There is no claim of unique local support loads, automatic heat allocation,
exact finite-step nonlinear energy/momentum conservation, automatic event
localization or cross-tree hydrodynamic coupling. The multi-tree entrypoint below
extends sustained contact to independent articulated trees.
This single-tree entrypoint admits at most 21 fixed-frame point contacts on one
existing admitted articulated tree.

G0/G1/G3 regressions advance 20 consecutive returned-state sliding and sticking
steps, check rotated/translated conveyor work and support release, preserve
redundant-support motion under contact permutation, and exercise resolved
impact, stopping, curved/interior penetration, cancellation and internal-couple
refusals. An independent supported prismatic-joint oracle verifies complete
implicit carrier recoil with and without fluid inertia and balances kinetic,
spring, gravity, damping and actuator work.

### Simultaneous sustained contact between articulated systems

`partitioned::articulated::contact::systems::step_with_contacts` advances
independently moving trees with one shared coupled midpoint iteration. The load
callback receives all `ArticulatedStage` states and generalized kinematics at
the same stage time, permitting cross-tree forces without stale partner states.
Each tree retains its own physical tolerances, dry-plus-fluid mass, analytic
fluid bias, implicit joint impedance and exact memory transition. The common
iteration budget is the minimum of the supplied tree budgets. The original
single-tree smooth, stop and contact APIs reuse this midpoint owner.

Contact rows retain only their incident tree blocks, combining both endpoint
wrenches first when a pair belongs to one tree. Cached inverse-mass responses
assemble the simultaneous half-step contact mobility; one final solve per
affected tree applies its accumulated generalized impulse. There is no dense
global mass assembly or sequence of pairwise contact projections. Full carrier
recoil and off-diagonal joint inertia remain active. Separate tree masses do
not imply hydrodynamic cross-tree inertia blocks.

The owner admits one to 128 trees, at most 4096 total links, 256 links/eight
joints per tree, and 21 contacts. Contact identities are canonicalized without
changing caller system indices. Every candidate is accepted or refused
together, preserving the aligned reaction, actual law/work and determinacy
checks above. Internal linear and angular impulses cancel across all incident
trees; prescribed supports contribute the reported external spatial impulse.

The mandatory `ContactSystemsPath` callback admits every configuration path,
continued witnesses, frozen contact frames and omitted-feature coverage before
publication. Actual nonnegative gaps and arithmetic-only incoming-speed
allowances remain required. Automatic CCD, friction-mode event localization,
unique redundant local reactions and exact nonlinear finite-step conservation
are outside this bounded entrypoint.

Regressions cover repeated two-tree support, a three-tree stack, system/contact
permutation, sliding transfer, a synchronized cross-tree spring, and refusal
atomicity. A separate analytic prismatic-joint oracle checks dry and fluid
recoil, spring/damper elimination and work under both system orders. A
16-step accelerating-current trajectory connects this owner to uniform-flow
added inertia, balancing affine fluid momentum and relative-fluid storage
against ambient source work.

### Mechanical travel stops with coupled recoil and support

`partitioned::articulated::stops::step_with_stops` advances the original coupled
midpoint state through bounded mechanical-stop events. Explicit lower/upper
travel coordinates, restitution and activation resolutions belong to each
joint. The shared normal-impact reset receives signed joint-coordinate rows
and the complete physical dry-plus-fluid mass, so simultaneous impacts retain
carrier recoil and all other joint responses. Impacts leave configuration and
aerodynamic memory unchanged; their mechanical loss is reported separately.

The original ballistic query locates closing and departing-return roots of
each candidate's quadratic joint path. Shortened trials are reintegrated from
the unchanged state. Only specifically typed native position-limit crossings
may shorten a refused trial; native limits are never removed or widened. A
contact state remains on the admissible side within its declared activation
band. No joint position is snapped to a limit. Smooth-step, minimum-interval,
attempt and event budgets bound the work and unresolved chatter refuses.

Exactly stationary active stops use the original SPD orthant solver on
midpoint velocities and the same effective mass as the smooth solve, including
its implicit spring/damper diagonal. Compressive reactions hold the joint;
noncompressive reactions release it. Other coordinates retain the full coupled
response. The reaction receipt corresponds to the velocity increment actually
accepted by the midpoint owner. Only arithmetic-roundoff residuals may become
exact zero rates, independently of user tolerances; impact cleanup is checked
again against complete generalized momentum and energy. All accumulated
reaction, work and loss outputs must remain finite. Memory advances only with
accepted smooth time, and a late refusal publishes no partial window.

`articulated_stops` checks independent Schur-inertia rebound and momentum,
simultaneous two-joint locking, stiff hold/release/recontact, a native-admissible
wide trial that departs and returns to a stop, and forced simultaneous support
with one joint releasing while the other retains its coupled reaction. Exact
memory timing and late atomic refusals are included. These are numerical
mechanical-stop models: `maximum_step_s` must resolve nonlinear smooth motion
and forcing. They are not interval-certified CCD of the exact trajectory,
stop compliance/damage, friction or historical control-mechanics validation.

## Hydrostatic pressure and material atmosphere loads

`atmosphere::MaterialAirProfile` composes the L2 atmospheric wind with the
existing L3 material gas owner at L4. It either retains one coherent arbitrary
gas snapshot or uses `StandardAtmosphere1976` with an explicit geometric ground
elevation. Wind heights remain above-ground distances; their sum with ground
elevation is the profile's geometric sea-level altitude. Queries preserve the
source's 0..80 km domain and retain its complete pressure, temperature, density,
sound speed and transport state. No duplicate EOS or hydrostatic law is added.

Frozen and evolving samples use the original mean/turbulence flow. They
transform both vector and gradient into flight-line forward/right/up, without
advancing an OU clock or consuming randomness. Legacy `AirScenario` constants
cannot override the material state. Homogeneous gas and the existing SimLoop
initialization/checkpoint behavior remain unchanged.

`OpenLoopDesign::force_buildup_in_atmosphere` and its evolving counterpart
subtract the sampled local wind from explicit inertial ground velocity, use
the shared longitudinal frame transform, and send this sample's density into
the original wing, drag and propeller load solve. G0/G3 tests compare the full
native result, verify hydrostatic `dp/dz = -rho*g`, pressure/density scaling of
dynamic pressure and Reynolds number, real lift/thrust scaling with altitude,
and exact read-only wind sampling. No extra constant reference density appears.

This is prescribed CG-sampled wind and a gas state frozen across the aircraft
for one load evaluation. The standard atmosphere is a dry-air reference, not
local-weather evidence. The original reduced longitudinal loads still omit
sideslip and Reynolds/Mach-dependent polar corrections; their gravity/frame
convention is unchanged. The bridge does not model pressure work, buoyancy,
weather evolution or validate high-altitude aircraft performance, and does not
make the historical SimLoop stratified. The prescribed solenoidal velocity
field does not imply divergence-free mass flux in a stratified gas; this is
not an anelastic or compressible atmospheric-flow solve.

## Six-DOF rotor momentum and gyroscopic motion

`spine::rotors::{RotorAssembly, RotorState, RotorTorques, step_held}` couples
up to 64 independent axisymmetric rotors to the native six-DOF spine.
`reduced_body` includes the whole assembly's mass and transverse/offset rotor
inertia but excludes rotor axial spin inertia. Its full tensor and all rotor
axes use the same declared body frame, which need not be principal. Signed absolute axial momenta
are `p = J (relative_spin + axis dot omega)`; `state_from_relative_spin`
constructs them from physical relative speeds. This is the reaction-wheel
reduction documented by Basilisk's `mrpSteering`, equations (2)--(6).

Motor torque changes rotor momentum and produces equal/opposite carrier recoil
without becoming external torque. Rotor aerodynamic torque changes both rotor
momentum and total external torque. `loads.moment_nm` excludes those rotor
aerodynamic torques, which the owner adds once. The native `fs-time` gyrostat
map includes precession from all stored rotor momenta. Translation uses a
world-force drift/kick/drift update, preserving its discrete orbital impulse.
The endpoint state includes unwrapped relative rotor phases and absolute
momenta; the receipt returns relative speeds, fixed-origin linear/angular
impulse residuals, separate motor/external work, actual kinetic energies and
the numerical work defect. No numerical defect is relabeled as dissipation.

The successful path is allocation-free, O(rotors) plus at most five elementary
gyrostat rotations, without a nonlinear solve. State, inertia, axes, loads,
count and finite derived quantities are admitted before publication. G0/G1/G3
tests in `tests/rotor_spine.rs` cover analytic motor recoil/phase/work,
external aerodynamic braking, multiple-axis impulse balance, balanced
counter-rotation and zero-rotor reduction, independent locked-inertia energy,
gyroscopic response, second-order phase/work refinement, replay and refusal.
The existing spine battery stays unchanged.

This is second order for the declared held loads, fixed axes and constant
inertia; coupled state-dependent loads need appropriate stage evaluation.
Physical rotor inertia must not be inferred from reflected engine/gearbox
shaft inertia. Gears, gimbals, imbalance, flexible rotors, changing mass and
contact are outside this owner. No historical calibration, complete rotor
checkpoint wire format or integration into the reduced longitudinal SimLoop
is claimed by this additional six-DOF API.

### Native passive shaft springs and bearing drag

`RotorAssembly::shaft_motion` prepares `passive::ShaftMotion` from the original
Flyer state and fixed `ShaftProperties`. It advances the same six-DOF carrier
through `fs-time::lie::gyrostat::rotors::TorsionalRotorStepper`, retaining the
native canonical momentum across accepted steps. Existing axes, absolute
shaft momenta and phases cross admission without normalization or a relative
speed round trip. The common inventory bound is 64.

Spring recoil, relative phase, spring storage and viscous heat come from the
native coupled solver. Zero stiffness selects its viscous path; zero stiffness
and damping select its held-torque path. Translation uses the original
drift/kick/drift map. Continuous force increments retain resolved addition
rounding in momentum and work residuals; nonfinite, underflowed and fully erased
transfers refuse. Explicit event impulses retain their original strict
transfer-relative admission. Motor work uses relative shaft motion, while
aerodynamic work includes absolute shaft motion exactly once. Physical heat,
elastic storage and the native gyrostat drift defect remain separate in the
receipt.

The prepared owner publishes its canonical and Flyer states together after all
motion and accounting checks succeed. Unresolved carrier momentum, phase or
impulse, invalid coefficients and native timestep-resolution failures refuse
without changing the accepted state. This canonical lane has its own precision
domain; the existing prescribed held-torque and instantaneous-impulse APIs
retain theirs. The analytical drag/oscillator, full-tensor trajectory, source
work and late-refusal regressions live in `tests/rotor_spine/passive.rs`.

`ShaftMotion::advance_with` evaluates a smooth carrier/motor/aerodynamic load law
on the accepted state and a native half-interval prediction of every carrier and
shaft coordinate. Both stages start from retained canonical momentum; the
midpoint load advances the original state through the full interval. Full- and
half-interval native factors are prepared once, so successful stepping adds no
allocation beyond caller callbacks. Only the accepted full interval supplies
work, spring storage and viscous heat; predictor transfers are never committed
or credited. This explicit midpoint composition retains second order for smooth
admitted laws, with refinement required for load feedback; it is not an implicit
coupling or event locator. Clock collapse, callback refusal and late accounting
failure leave both stored states unchanged; callback side effects remain the
caller's responsibility. `tests/rotor_spine/passive/stages.rs` checks native
full-tensor stages after impulses, independent work/momentum accounting, a
coupled ODE refinement oracle, exact held-law parity and refusal/retry.

### Real aerodynamics through independent passive shafts

`propulsion::passive::ShaftPropulsion` binds the existing two-disk wing/BEMT and
native lifting-surface owner to `ShaftMotion`. Each mounted shaft supplies its
own signed absolute speed `p/J`; blade handedness maps that speed into the
existing positive-spin BEMT domain. Axial aerodynamic torque enters only the
external shaft port. Motor torque remains internal, while force-arm and
transverse hub moments plus world gravity use the original carrier ports.

Axial and spatial disk-wash paths, including custom lifting surfaces, share the
same aerodynamic stage as the geared and dynamic consumers. Two actual solves
sample the initial and native predicted midpoint states, including after
impulses; only the accepted full interval contributes mechanical work and heat.
A physical-inventory mismatch refuses before callbacks. Any stage or final
accounting refusal preserves both native and published motion; callback side
effects remain caller-owned.

This consumer is quasi-steady. It adds no dynamic wake storage, ideal gear
constraint or total fluid/body energy closure claim. `tests/passive_propulsion.rs`
checks axial/canted shared-owner parity, independent rotor inventory mapping,
real two-stage stepping after an impulse, independently reconstructed physical
momentum and work, refusal/retry, and the distinct numerical admission domains
of continuous forces and explicit event impulses.

`propulsion::passive::SwirlShaftPropulsion` explicitly selects the existing
native annular axial/angular momentum law on the same independent shaft owner.
It shares the geared swirl geometry dispatch and the passive admission,
force/torque assembly and transactional midpoint advance. Each shaft's actual
signed absolute rate drives its own accepted rotor solve; aerodynamic torque
enters that shaft's external port exactly once, including counter-rotation and
regenerative windmilling. The accepted annular solutions remain attached to
the same forces and torques in `ShaftEvaluation<N, CoupledSwirlStep>`.
The default payload of `ShaftEvaluation` and `ShaftPropulsionStep` remains
the ordinary axial `CoupledStep`.

Only a complete accepted aerodynamic stage supplies the next numerical warm
seed; it adds no physical wake inertia. Springs, bearing heat, motor work and
absolute aerodynamic shaft work retain the native mechanical owner. The G0/G1/G3
regressions in `tests/passive_propulsion.rs` target identical-state axial/canted
geared parity, independent post-impulse midpoint momentum and energy accounting,
zero-motor windmilling into free shaft storage, bounded solver refusal and
midpoint rollback/retry. Wing and mutual-disk wash remain the original axial
actuator-disk field; no dynamic swirl, helical wake transport, or whole-fluid
energy closure is added.

### Dynamic wakes through independent compliant shafts

`propulsion::passive::dynamic::DynamicShaftPropulsion` composes the same prepared
`ShaftMotion` with the existing constant-density uniform axial dynamic-inflow
owner. One private `DynamicWake` owns the original wing-at-induction, blade,
mutual-wash and implicit wake-flux calculations for geared and independent-shaft
consumers; existing harmonic consumers route through that same implementation.
Each mounted propeller receives its absolute `p/J`. Internal motor and external
shaft torque, carrier wrench, spring storage, bearing heat and retained
canonical momentum keep their original mechanical owners.

`DynamicShaftState` supplies controls with the actual joint stage. Initial loads
predict mechanics and both wakes; midpoint loads advance the original full
interval. Both wake trials finish inside the mechanical load callbacks, so an
initial, midpoint, full-wake or final-mechanical refusal leaves canonical motion,
published state and supplied induction unchanged. Successful publication ends
with an infallible induction assignment; callback side effects remain caller-owned.
Optional oblique quadrature and pointwise material-air modes reuse the existing
native implementations.

Apparent air mass is not aircraft mass. Per-disk axial impulse, storage and work
are reduced-state accounts, not resolved-wake world momentum or whole-fluid/body
passivity. Independent-shaft harmonic/skew wake storage and unsteady wing memory
remain outside this consumer. Six cases in `tests/passive_propulsion/dynamic.rs`
check shared-owner parity, post-impulse native stage oracles, independent physical
source accounts, joint second-order refinement, shaft-axis invariance and late
refusal/retry.

### First-harmonic and skew wakes through independent shafts

`propulsion::passive::dynamic::first_harmonic::HarmonicShaftPropulsion` adds
the original physical mean and two first-harmonic modes per disk to retained
`ShaftMotion`. A private `HarmonicWake` shares prepared geometry, pointwise
blade/wing evaluation, modal admission, normal/skew transport and coupling-work
accounting with the existing geared and continued harmonic consumers. Native
modal laws and their bounded skew solver are unchanged. Each independent rotor
still supplies its actual signed `p/J`; no ideal gear is introduced.

The six induction velocities retain their fixed body disk frames through
changing incident wind. Carrier frame-rate terms remain distinct from physical
shaft speed. Oblique hub forces/moments, axial aerodynamic shaft torque, apparent
fluid storage, signed skew cross work, spring energy, bearing heat and motor
work keep the original owners and enter each accepted account once.

Initial loads predict both canonical mechanics and every wake mode. The full
modal/skew proposal and finite coupling-work defect are computed inside the
second load callback, before mechanical publication. Any aerodynamic, wake,
coupling-accounting or late native mechanical failure preserves both stored
mechanical states and all caller-owned modes; acceptance ends with an infallible
mode assignment. Predictor transfers are not accepted work or heat.

Normal/skew snapshot parity, asymmetric post-impulse stage oracles, independent
mechanical and modal momentum/work, signed transport decomposition and late
refusal/retry are checked in
`tests/passive_propulsion/dynamic/first_harmonic.rs`. Smooth coupled trajectories
still require timestep refinement; no whole-fluid momentum/passivity theorem,
dynamic swirl, transported wake geometry or unsteady wing memory is claimed.

### Continued wing branches through independent shafts

`HarmonicShaftPropulsion::initialize_continuation`, `evaluate_continued` and
`step_continued` reuse the original `fs_wing::nonlinear::continuation` owner
with the same normal/skew modal wakes and retained `ShaftMotion`. The caller
keeps the native `WingContinuation` alongside the physical modes; it stores an
accepted algebraic circulation/incidence branch, not added fluid inertia or a
lagged aerodynamic force. Plain and continued paths share one stage evaluator
and one joint mechanical/wake stepping implementation.

The initial force solve uses the caller's accepted anchor, and the midpoint
solve uses the newly solved initial anchor. Both stages recompute the actual
wing/blade loads at current pose, independent signed spins and all fixed-frame
induction modes. The accepted midpoint anchor becomes the next initial solve's
guess. Native geometry identity, incidence continuity, row conditioning,
filament clearance, residual and bounded solver-work gates remain unchanged.

All proposed anchors are local until the full step succeeds. Each stage's
returned anchor is validated before its load callback returns; the full modal/
skew wake and coupling-work checks also finish before canonical mechanics can
publish. Any branch, material-air, wake, accounting or native mechanical refusal
preserves canonical and published motion, all six modes and the accepted wing
anchor. Success moves the already admitted modes/anchor into place without a
new fallible operation; callback side effects remain caller-owned.

The focused cases in
`tests/passive_propulsion/dynamic/first_harmonic/continuation.rs` exercise
48 ms normal and skew trajectories after shaft impulses, exact stage/anchor
replay, independently reconstructed physical impulse and work, modal storage
accounts and joint timestep refinement at 6/12/24/48 steps. Body, independent
shaft and modal endpoint differences, together with coupling-work defects,
show second-order refinement within that declared fixture. Late branch,
material-air and final native-mechanical refusals preserve every accepted state
and replay exactly. This extends algebraic branch continuation; it adds no
unsteady-wing constitutive law, free-wake geometry or whole-fluid passivity
claim.

### Ideal gears and real propeller/airframe motion

`spine::rotors::geared` admits one ideal relative gear coordinate and 1–64
axisymmetric rotors with explicit signed ratios `g_a`, physical axial inertias
`J_a`, and body-fixed unit axes `e_a`. The caller supplies the whole assembly's
mass and full locked inertia tensor, including all rotor inertia. Assembly
computes `J_eq = sum(J_a g_a²)`, `c = sum(J_a g_a e_a)`, and the SPD Schur
reduction `I_red = I_locked - c cᵀ/J_eq`. Products of inertia are retained by
`fs_time::lie::gyrostat::InertiaTensor3`; its principal-axis assembly occurs once.
State is common relative phase `s`, canonical momentum
`p = J_eq s_dot + c.dot(omega)`, and the full carrier six-DOF state.
Total body momentum is `I_red omega + c p/J_eq`. Physical kinetic energy is
`m |V|²/2 + omega.dot(I_red omega)/2 + p²/(2 J_eq)`; the battery independently
reconstructs it from the original locked inertia and relative rotor rates.

Internal generalized motor torque changes `p`; aerodynamic torque changes it
by `sum(g_a tau_a)` and contributes `sum(e_a tau_a)` to total external torque
exactly once. Carrier moments must exclude those axial shaft torques. Motor
work is `Q_motor delta_s`; aerodynamic shaft work uses absolute axial rate
`g_a s_dot + e_a.dot(omega)`. A perfectly balanced counterrotating pair has zero
`c` and no motor recoil; unequal external shaft torque still changes total
angular momentum. Signed ratios are never replaced by reflected inertia in
that cancellation law.

`geared::step_held` is second order for fixed world force/body torques;
`geared::step` evaluates a smooth load law at the initial state and a half-step
predictor before accepting a full step with the midpoint loads. The two calls
occur in deterministic start/midpoint order. The work residual is the actual
endpoint kinetic change minus declared work, with no artificial dissipation.
The global angular impulse combines the translation kick's orbital impulse
with trapezoidal endpoint-frame torque impulses; its residual measures the
accepted discrete momentum balance. Successful mechanics use fixed-size
storage. Invalid models/states refuse; input state is immutable, while side
effects inside caller-owned callbacks remain the caller's responsibility.

`spine::rotors::propulsion::GearedPropulsion` is a concrete consumer: both
stages run the existing distributed two-way wing/BEMT solve with actual stage
attitude, velocity, all three angular rates and generalized shaft momentum.
World-position material-wind samples are converted to body relative incident
velocities at every strip and both disk centers. The converged circulations
and accepted local wash drive full vector Kutta–Joukowsky forces, profile drag,
section moments, existing classical induced drag and parasite ledger, and
thrust at both disk arms. Gravity is added once in world axes. Force moments
retain all three components about the same assembly COM. Propeller mounts can
refer to two rotors inside a larger gear train, so additional engine/flywheel
inertia is supplied physically by the caller rather than invented.

The shared steady, dynamic, harmonic and swirl load assembly now streams its
accepted local strip velocities and panel circulation through
`fs_wing::wrench::CirculatoryWrenchAccumulator`, with the assembly COM as the
moment origin. This returns the complete circulatory force and moment through
the common spatial-wrench type. The original strip/panel order and Vec3
operation order are retained bit for bit; profile drag, section moments,
induced/parasite drag, thrust and shaft-torque ownership remain separate.
Invalid or overflowing panel wrenches refuse before a load evaluation is
published. The new local-flow hinge API uses that same panel force/moment law.

The native loop migration is covered by exact six-component arithmetic parity
and independent spatial-work/origin/covariance tests. The 29-test actual wing,
hinge and nonlinear numerical-source run passes, including the existing hinge
golden; the complete native propulsion and frontend suites were not executed
for this change. No historical SimLoop path or solver equation is changed.

Gear ratio, rotor axis, blade handedness and positive thrust direction are
distinct inputs. Handedness maps signed ABSOLUTE spin into the current BEMT
positive-speed branch; the physical aerodynamic torque reverses with that
handedness. The current actuator-disk coupling only supports body-x shafts and
explicit +x thrust. Stopped/opposite blade rotation returns
`propeller-spin-outside-domain`; reversed net disk flow retains the existing
`disk-inflow-invalid` refusal. Neither spin nor inflow is floored. Assembly
and stage defects use `propulsion-assembly-invalid`, `propulsion-state-invalid`
or existing native refusals. No unconverged/one-way aerodynamic fallback runs.

G0/G1/G3 evidence in `tests/geared_propulsion.rs` covers nonunit gear-ratio
recoil/phase/work, balanced and unequal torque cases, arbitrary-axis dense
inertia and independent physical momentum reconstruction, admitted-count and
state boundaries, current symmetric aerodynamic parity, real asymmetric
six-DOF motion, absolute-speed/handedness refusal, and second-order refinement
of both an analytic motor/drag law and real coupled aerodynamic dynamics.
The `GearedPropulsion` aerodynamic tier remains quasi-steady, with existing 2-D section closures,
axial disk inflow, straight wakes, spanwise wash overlap and Estimated drag.
It does not add rudder dynamics, propeller crossflow/harmonics, apparent mass,
shaft compliance/failure, gimbals, contact, checkpoint serialization, or
historical calibration. Historical preset mass/inertia and SimLoop behavior
are unchanged. Full coupled accuracy remains bounded by the declared wing,
BEMT and coupling tolerances; refinement evidence does not certify physical
validity or arbitrary operating branches.

### Wind-axis lift in aircraft loads and native propulsion

The longitudinal aircraft and native steady, dynamic and harmonic propulsion
consumers share one projection helper to compute the
circulatory lift component along `(sin(alpha), 0, -cos(alpha))` before applying
the existing classical induced-drag screen. Body-normal force `-Fz` is only its
projection; squaring that component would spuriously reduce induced drag by
`cos(alpha)^2` in uniform flow. The correction changes neither the section
solve nor any profile-drag, parasite, propeller or moment-arm owner.

The local analytic regression specifies signed-lift recovery under pitch,
streamwise-force and side-force changes. The real geared-aircraft tests
independently integrate the accepted flat-panel circulation and use a velocity
projection to compare the corrected native force at attached and separated
states, retaining the existing force/moment tolerances and shaft-power checks.
This remains an Estimated global induced-drag model, not a resolved-wake drag
calculation. `ForceBuildup::lift_n` retains its documented circulatory body-up
load; the transverse force is used for induced drag before profile/parasite
forces are added. The separated aircraft regression independently reconstructs
lift from accepted local flow, circulation and span, verifies induced drag and
actual drag work, and exposes the body-axis error at 0.5 rad incidence. In
uniform flow that deficit is 23%. The longitudinal SimLoop
and trim now receive this same corrected force; their trajectories and physical
trim fingerprints intentionally change, without changing the section/coupling
solve or the model's longitudinal applicability.

### Instantaneous rotor and geared momentum transfer

`spine::rotors::apply_impulse` and `geared::apply_impulse` apply simultaneous
carrier, internal motor and external axial rotor impulses without advancing
time. `CarrierImpulse::at_body_point` forms the COM angular impulse from a
body-frame arm and world-frame linear impulse. Full-tensor carrier recoil and
signed gear ratios preserve total world angular momentum under internal
drive; external shaft impulses enter the total exactly once. Motor work uses
mean relative rates, while external axial work uses mean absolute rates.

Positions, orientation and shaft phases stay fixed. Dynamic and harmonic
consumers replace only their public mechanical `motion` state, retaining wake
induction and histories until the next actual force evaluation. No contact,
restitution or fluid-added-mass impact law is inferred. Per-coordinate transfer
checks refuse `impulse-unrepresentable` when the state cannot resolve the jump;
large background momenta cannot hide a dropped impulse. Increment-based
momentum and polarized kinetic-work balances retain small physical changes
without subtracting large endpoint totals. The `rotor_spine` battery covers
analytic recoil, signed gears, off-center impulses, world-frame covariance and
representability with a large neighboring rotor momentum.

Finite held-load updates also refuse a completely erased nonzero shaft kick,
torque-times-step underflow, or a shaft change lost when forming the aggregate
stored rotor momentum used by the gyrostat. The shared geared midpoint owner
carries those refusals into native steady/dynamic/harmonic propulsion. Normal
finite-step addition quantization near balanced motor/aerodynamic torques is
retained; the stricter transfer-relative budget applies only to an explicit
instantaneous jump. Accepted finite trajectory arithmetic is unchanged.

### Spatially mounted propellers

`propcoupling::spatial::{DiskAxes, coupled_step}` and
`GearedPropulsion::new_spatial` extend that consumer to canted, body-fixed
shafts. Positive-thrust axes are admitted unit vectors, collinear with their
physical rotor axes; neither gear sign nor blade handedness chooses the
thrust direction. Native wing-induced velocity and the actual local rigid-body
incident velocity are projected onto each disk's axis before BEMT. The returned
thrust acts along that axis at its actual arm, while aerodynamic shaft torque
acts along the separately signed spin axis, exactly once in total momentum.

Spatial wash measures radius in the full disk plane, including vertical
separation. Each disk's scalar induced velocity becomes a vector increment
along its thrust axis at every strip; circulation and full-vector force
assembly use the same accepted local velocities. Geometry factors are prepared
once per coupling solve. The existing Aitken solver, per-disk convergence gate,
momentum root and refusal behavior are shared. The spatial recipe and actual
axes have a distinct content digest; `new` and the historical axial APIs keep
their original spanwise model and numerical behavior.

G0/G1/G3 checks in `tests/spatial_propulsion.rs` independently reconstruct
disk-plane overlap, projected native BEMT thrust/torque, vector wash and global
angular momentum; verify rotation covariance, the coplanar axial limit, and
invariance to reversing a mechanical axis and its coordinate signs; and reject
invalid mounts or a late failed stage without returning partial motion. This
remains a disk-average, quasi-steady axial model. Transverse disk loads,
azimuthal harmonics, swirl, contracted/skewed wakes, direct disk-to-disk wake
interference and gimbal motion are absent. The finite 2R wash window remains
an explicit model truncation, not a resolved far wake. Canted-mount support
does not establish a crossflow validity envelope or historical calibration.

### Explicit native lifting surfaces and section frames

`propulsion::surfaces::LiftingSurfaces::new` admits borrowed panels, strips,
per-strip `SectionProperties`, and declared `InducedDragGroup` references.
`GearedPropulsion::new_spatial(...).with_lifting_surfaces(&surfaces)` selects
that body-fixed layout in every steady, native-swirl, dynamic, harmonic and
explicitly continued aerodynamic stage. The shared material-air sampler,
physical wake transport, rotor/shaft ports and midpoint mechanical owner are
retained. The original design still supplies rotors, disk positions, CG and
the parasite ledger. Selecting no custom layout preserves the built-in path.

Each admitted panel has exactly one strip owner; every strip has one declared
drag group, and every group must be used. Bound rows in a strip share their
physical plane, directed span, width and surface identity, with chordwise
registration along the actual panel chord, including swept streamwise rows.
Frames and true bound lengths are prepared once. The supplied strip chord is
the physical chord normal to the bound span. Base strip twist receives
`canard_gain * canard_rad` once; this is a section-polar incidence control,
not physical panel or hinge articulation. Fixed incidences and zero control
borrow the original strip layout without cloning its panel-index vectors.

Complete local air and vector propeller wash reach the framed nonlinear owner
before section projection. Full velocity and actual bound segments supply
Kutta–Joukowsky force. Profile drag uses the accepted effective incidence,
section-plane speed/direction and chord times true bound length; spanwise
velocity cannot inflate the two-dimensional dynamic pressure. The camber
pitching free couple follows the local span, including yaw couples for a
vertical surface. Custom layouts require finite CG speed above 1 m/s and
native local chordwise/disk-flow admission, without a hidden global body-x
velocity or body-xz incidence restriction.

Induced drag remains an explicitly declared classical engineering closure.
For each group, sum only its circulatory forces, remove the component along
CG relative wind, and apply `D = |F_perp|² / (q_CG * pi * e * span²)` opposite
that wind at the CG. Reference area cancels algebraically. Grouping related
biplane surfaces under one effective reference is distinct from assigning
separate references to independent or opposing loads. No cross-group
induced-drag interference, extra drag moment, or Trefftz-plane accuracy is
claimed; mutual induced velocities remain in the coupled circulation solve.

`tests/native_surface_wrench.rs` independently reconstructs complete canted,
swept and vertical loads, grouped drag and a rotated local-forward scene.
`tests/section_surface_stages.rs` checks the accepted section equation in all
native wake lanes, explicit incidence-control ownership, and actual midpoint
orbital/carrier/rotor momentum balance with vertical-surface forces.

## Canard hinge friction and stops

`canardmech::step` solves one monotone scalar midpoint-velocity equation
with a bounded safeguarded Newton/bisection iteration. A successful step
uses that same velocity for friction, stop damping, and applied work;
`StepReceipt` exposes all three. The quadratic stop spring uses its discrete
energy gradient across engagement/release, and stop damping acts over the
fraction of the configuration path beyond either stop. Mechanical energy
change plus both losses equals held-drive work to solve/roundoff tolerance.
The smooth Coulomb law still permits sub-Coulomb creep; exact static sticking
is not claimed. Non-finite derived steps and unresolved bounded solves refuse
with `mech-step-non-finite` and `mech-step-unresolved`, respectively.

## Invariants

- Admission precedes every computation; an unadmitted design cannot produce
  a build-up or panel.
- The declared empty mass must equal the component sum within 1 kg
  (`mass-spec-mismatch` otherwise — no hidden mass, ever).
- The hinge-axis admission domain IS the E1.5 prior [0.25, 0.50] x/c.
- Reference-config positions/gyrations are a CALIBRATED reconstruction
  (documented in-source): positions place the gross CG at the dossier's
  29.7% chord; gyrations land the inertias inside the published ±15% band.

## Spatial wind in the flight lifecycle

`SpatialWindSpec::spectrum` explicitly selects `Analytic { sigma_mps,
length_m }` or `MannStress { parameters, reference_height_m }`. The latter
fits the existing wall-compatible basis to the target's u²/v²/w²/uw at one
height, then evolves the same paired OU quadratures used by every spatial
constructor. The calibrated velocity reaches the CG or actual wing/canard
strips and propeller disks, including distributed roll/yaw loads. Fitting and
Mann quadrature occur at initialization; per-tick sampling uses the original
allocation-free modal evaluation. `spatial_wind_stress_fit()` exposes the
achieved finite-basis residuals. This is one-point stress calibration, not a
complete spectral/coherence fit or a recovered historical gust trace.

Calibrated runs use V7 checkpoints, retaining the target parameters, reference
height, mode count, time scales, mean, advection, both OU states, sampling
mode and optional material gas/profile. Restoration rebuilds the exact
seed-specific fit and checks run identity before accepting the state. The
original analytic model retains its V3–V6 bytes and realization. G0/G3/G5
tests in `simloop/spatial_wind/mann_tests.rs` cover the exact atmosphere-to-
flight field, actual rotor/ground/lateral response, all six sampling/material
combinations, model-substitution refusal and self-contained replay.

`MannHomogeneousAbove { parameters, min_kl, max_kl, minimum_height_m }`
connects the complete finite-band homogeneous Mann field to the same native
CG, wing/canard strip, propeller and distributed lateral load owners. Every
wind query requires height at least the positive authored minimum. This is an
applicability restriction for exploration, **not proof that ground effects
vanish above that height**; the model has no wall projection or surface blocking.
Its mode count comes from the explicit correlation-time vector, and omitted
wavenumber bands are not renormalized. `homogeneous_wind_covariance` reports
the represented turbulence-only space-time tensor in flight coordinates.

Homogeneous runs use V8 checkpoints with target, cutoffs, applicability height,
all four OU quadratures, material choice and sampling mode. Tick-zero identity
contains every quadrature; restoration rebuilds and checks the same model.
The V3–V7 paths retain their bytes and behavior. A domain refusal at a CG,
strip, disk, midpoint or accepted endpoint returns without consuming flight,
pilot/perception or random-path state. An opt-in rollback buffer is allocated
once and reuses control rings and OU arrays on subsequent steps; immutable
spectral factors are not rebuilt or copied each tick. The homogeneous replay,
load, covariance, model-substitution and atomic-refusal regressions are in
`simloop/spatial_wind/homogeneous_tests.rs`.

`MannPlaneBlocked { parameters, min_kl, max_kl }` applies fs-atmo's existing
harmonic plane correction to the same full-tensor Mann field at the aircraft's
ground h = 0. Every CG, wing/canard strip and propeller sample uses that corrected
velocity and analytic gradient. Normal wind vanishes at the ground; the field
stays solenoidal, retains all original modal cross-components, and approaches
the homogeneous field aloft. The mean log law and material gas/profile continue
to use their original owners. This is an explicitly selected outer-flow slip
model, without a viscous sublayer, ground-generated turbulence, a near-ground
spectral calibration or reconstructed historical gusts.

`mann_wind_covariance(first, second, lag)` exposes the complete represented
space-time tensor for either full-tensor model in flight coordinates. Plane
blocking acts at both observation points, so equal separations at different
heights can have different covariance. Sampling retains the same four sequential
OU paths; neither sampling nor covariance draws noise or advances the wind clock.
Plane-blocked runs use V9 checkpoints with a distinct model identity, target,
cutoffs, four paths, sampling mode and material choice. V2–V8 keep their existing
layout. Below-ground material-point refusals reuse the existing transactional
flight/control/OU backup. Focused integration regressions in
`simloop/spatial_wind/plane_blocked_tests.rs` cover wall parity, all field/tensor
components, actual load and motion changes, replay and refusal rollback.

Ordinary, CG spatial-wind and distributed-wind airborne paths share one coupled
longitudinal/rotor midpoint step. The first native aerodynamic solve predicts
midpoint velocity, pitch attitude/rate and rotor speed. A second solve evaluates
the actual loads at that predicted state, including spatial wind at the stage's
CG and material points. Both stages use the same held controls and temporal wind
realization; the legacy scalar gust's incidence increment uses stage airspeed.
The first solution warm-starts the second, and only the accepted midpoint wash
is cached for the next tick. There is no duplicate aerodynamic owner or extra
atmosphere/control clock advance.

The force point uses `X_mid = X_n + dt V_n / 2`; accepted centre-of-mass drift
uses the mean of initial and final ground velocities. This closes both linear
impulse and orbital angular impulse about any translated origin to roundoff:
`Δ(X × mV) = dt X_mid × F_mid`. Pitch and rotor inertia receive their respective
midpoint torques. Smooth coupled longitudinal/rotor dynamics are second order
under the held inputs; this does not extend the claim to OU innovations,
discrete controls, existing load-envelope clamps or the reduced lateral model.
Constant force/gravity and constant angular accelerations integrate exactly.
A refused midpoint publishes no predicted motion or warm cache. A positive
midpoint CG reaches the native force solve even when the held-load endpoint
predictor would cross ground. Actual accepted midpoint kinematics determine
GroundContact; the terminal state retains that midpoint velocity, pitch/rotor
impulse and warm cache, with height placed on the wall. Distributed material
points remain subject to native atmosphere admission: a positive CG does not
certify that every surface is above ground. Only a nonpositive midpoint CG uses
the bounded held-load terminal fallback without requesting loads underground.
These event branches claim no smooth order or continuous collision/localization.

The ordinary path retains body airspeed as stored state, converting through the
steady headwind offset before and after the inertial step. Force-free attitude
changes therefore preserve world momentum and kinetic energy to roundoff instead
of creating energy through explicit rotating-frame transport. Airborne ground
motion may be backward. Both paths remove the force assembly's level-flight
gravity convention and apply vertical world gravity once, independent of angle
of attack, pitch attitude or a gust. Focused analytic regressions live in
`simloop/motion.rs`. Ordinary and spatial launches also share one constrained
world-momentum step: stationary pitched aircraft have both body airspeed
components, the loaded rail holds vertical ground velocity at zero, and actual
rotated non-gravity loads determine support, friction and release. If tick-start
upward load meets or exceeds weight, the same full airborne midpoint interval
begins immediately, including vertical/pitch/rotor and lateral evolution, with
no preceding rail friction or tow impulse. Release introduces no velocity jump
and advances controls/atmosphere only once. Positively loaded rail ticks retain
the existing single-force-solve constraint/tow law; within-tick loss of support
is not localized.

Pilot member 3 retains its registered gains. Its earlier crash-class calibration
predates these momentum/gravity corrections and is no longer a validation claim
for the current plant. With consistent stationary/rail velocity, the unchanged
seed-1903 Dec-17 reference reaches RailEndWithoutLift at tick 975, approximately
18.314 m: the 18.3 m rail is too short for this currently calibrated plant.
The battery pins that outcome rather than fabricating lift or retuning gains.
Airborne lifecycle tests explicitly add the existing registered catapult as a
test fixture, distinct from historical Dec-17. With consistent disk-induced
wash and same-tick support release, member 3 flies 11.42 s with five undulations
and 44.90 m of travel before the native wing solver exhausts its 600-evaluation
budget at tick 1444 (residual 5.92825e-10). Fixed controls fly 4.84 s, reach
11.61 m forward and turn back to 5.41 m before the outer coupling exhausts its
12-correction budget at tick 655 (residual 0.00584813, tolerance 0.005).
The native strip-cycle safeguard preserves the unchanged Huffman 4 m/s catapult
liftoff gate: lift occurs at tick 197 and 10.229 m, followed by 4.66 s airborne
before an outer-coupling second-growth refusal at tick 756. The unchanged
3 m/s case passes its initial wing solve but reaches the 120-evaluation
propeller station cap at tick 8. These are numerical nonconvergence receipts
carried by EnvelopeExceeded, not established physical envelope boundaries.
The substantive flight duration/undulation and travel gates remain, along with
exact terminal code/tick and replay evidence. Neither assisted fixture nor
original rail shortfall establishes a revised historical launch or pilot calibration.

`SimLoop::init_with_spatial_wind(scenario, SpatialWindSpec)` opts into the
wall-compatible evolving `fs-atmo` field. The specification declares the mean
log law, spatial strength and length, advection speed, and every modal OU
correlation time; the scenario supplies the seed and load density. Its mean
reference speed must match the scenario headwind. This finite-mode spectrum
is synthetic and is not a validated reconstruction of historical weather.

The longitudinal force and propeller solves consume CG-sampled body air-relative
velocity. Ground velocity is separate inertial state: temporal OU changes and
spatial wind gradients change relative airspeed without fabricating an aircraft
impulse. Forces use the accepted tick-n field; the returned tick-n+1 state and
read-only `sample_spatial_wind` queries share the advanced field. The old scalar
OU angle increment is disabled in this mode. Gravity acts vertically in world
coordinates at the actual pitch attitude. The sampler exposes the full velocity
and analytic gradient in flight-line-forward/right/up coordinates, but the
current longitudinal force interface consumes only forward and vertical inflow;
distributed span/disk gust loads and crosswind/sideslip dynamics are not claimed.
Mechanical initialization retains the baseline steady trim, not gust-equilibrated
prelaunch dynamics.

Spatial runs save V3 checkpoints containing the entire wind specification,
both actual OU quadratures, ground velocity, and run identity. The existing
`restore_checkpoint` resumes this state exactly and checks the relative-air
velocity against the restored ground/wind state. Legacy initialization and the V2
checkpoint layout retain their existing encoding. Focused tests live in
`simloop/spatial_wind/tests.rs`: actual load/propeller response, zero-amplitude
gust exclusion, observer independence, rail/airborne replay, and gravity-frame
consistency. Run with `cargo test -p fs-flyer --lib simloop::spatial_wind`.

### Distributed strip and propeller inflow

`SimLoop::init_with_distributed_spatial_wind` samples every wing/canard strip
reference point and both propeller disk centers from the same frozen tick.
`OpenLoopDesign::force_buildup_with_incident` admits the corresponding explicit
body-relative incident field. Ordinary `force_buildup` and `force_buildup_warm`
use this same owner with a uniform ambient field, including at zero pitch rate.
Rigid pitch velocity `q cross (position - CG)` enters each local sample once;
the former clamped pitch-as-twist approximation is removed. Both the axial
`q*(z-z_CG)` and vertical `-q*(x-x_CG)` components change actual local flow,
including the separate propeller-hub inflows. Local chord-plane speed and
incidence drive the nonlinear strip solve, and its accepted local velocity
drives Kutta-Joukowsky forces,
profile drag and section pitching moments. Converged propeller wash is included
in both circulation and force velocity. The separate disk-center axial inflows
drive the two BEMT thrust/torque solves and the existing engine dynamics.

Both native propeller/airframe paths carry the same disk-induced axial velocity
`vi`. Ideal actuator-disk momentum gives mass flux `rho A (V + vi)`, far-wake
increment `2 vi`, and thrust `2 rho A vi (V + vi)`. The rationalized disk root
therefore returns `vi`, not `2 vi`. Each strip receives `wash_factor * vi` once;
the factor already supplies the disk-plane value 1 and the upstream/downstream
spatial variation. The former additional factor 1/2 in ordinary loading was
inconsistent with its own momentum root. Both loading modes now use the same
wash sum and root. The coupling spec digest uses v2 to distinguish the corrected
physical convention; candidate tuples, acceptance tolerances and budgets are
unchanged. The independent native regression reconstructs mass-flux thrust and
the washed wing circulation, and rejects a second halving of the wash.

The compact uniform coupling entry also uses washed local incidence, chord-plane
speed and force, and returns `Some(strip_velocity_mps)` from its accepted
evaluation. These velocities reconstruct the returned circulation exactly; the
next `w_slip` predictor must not replace them. Uniform disk inflow is the actual
incident-plus-wing-induced axial velocity, with no 0.1 m/s floor. Small positive
inflow reaches native BEMT unchanged, zero remains admitted, and reversed or
nonfinite inflow refuses. Its result identity is scoped to
`uniform-prop-coupling.v1`; the already-correct explicit, spatial and framed
identities are preserved. The compact coupling golden changes to `0a50ecad...`
under the unchanged 350 RPM fixture. Independent mass-flux, accepted-flow,
uniform/explicit equivalence, 0.05 m/s BEMT and actual vortex-induced reversal
tests cover the correction. The scalar-density native 120-tick startup scenario
retains `8d87ebbb...`; material and WebAssembly qualification are separate.

The aircraft force assembly always evaluates Kutta–Joukowsky and profile loads
with those same accepted local strip velocities. Ordinary zero-rate loads and
trim also use this consistent local field, keeping the force law continuous
through zero pitch rate. `aircraft_battery` checks uniform/distributed identity
for signed pitch rates, assembled/solver lift agreement, finite-differenced
moving-point velocities, actual BEMT thrust/torque at both moving hubs, and
separated-drag work with strip-local pressure and force directions. The finite
2R influence window, spanwise overlap approximation and omission of wake contraction, swirl
and vertical disk overlap remain explicit model limitations. References: NASA
Glenn, [Propeller Thrust](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/propeller-thrust/),
and Deters, Ananda and Selig, AIAA 2015-2265, equations (3)–(7).

Both native propeller/airframe paths judge residual growth against the evaluated
predecessor of the current correction. A first growth above 25% rejects the
entire trial and retries `x_previous + omega_min * residual_previous`; the
rejected state never becomes the retry's origin or Aitken history. Subsequent
corrections retain the relaxation floor, and a second genuine growth refuses.
Both the rejected trial and retry count against the unchanged candidate cap;
there are at most `cap + 2` native map evaluations including initialization.
The same acceptance tolerance remains in force. Correct disk wash resolves the
former member-3 tick-1443 force input within the unchanged candidate-B tolerance;
that local recovery does not establish a complete flight envelope. The native
zero-cap fixture still refuses unresolved joint loads. No one-way fallback is
introduced.

Within each solve, wing evaluations now continue from the evaluated predecessor's
accepted circulation. The first wing evaluation retains its original cold seed;
subsequent accepted trials transfer the already-owned gamma buffer without a
copy. Rejected candidates cannot replace that anchor, so the retry restores both
the slip and its wing-branch history. This avoids cold-start basin changes being
mistaken for the physical coupled map. Continuous section-regime transitions are
allowed; no alternate-root search is added. The coupling spec identity is v3.
This policy does not retain wing circulation across ticks or prove uniqueness.

`propcoupling_continuation` freezes the former fixed-control tick-553 input.
Cold starts there alternated between valid wing roots producing roughly 4.361
and 3.209 m/s of disk induction. The continued branch closes after initialization
at a joint residual `7.216e-5`, below the original candidate-B `0.005` gate.
Independent trailing-field/section reconstruction gives a wing mismatch
`5.608e-11`, below its unchanged `1e-10` gate; native BEMT and actuator momentum
also close at the returned circulation and actual evaluated wash. The former
off-trim branch-chatter fixture likewise resolves within both original caps.
The actual fixed-control lifecycle now reaches native station nonconvergence
at tick 577; its 120-evaluation refusal remains explicit. The assisted member-3
flight still reaches physical ground contact at tick 1102 with three undulations.

`coupled_prop_airframe_step_with_incident` carries disk-induced velocity `vi`,
with `T = 2 rho A vi (V_disk + vi)`, and applies `wash_factor * vi` at each
strip. Returned strip velocities belong to the actual evaluated iterate that
produced circulation. Distributed acceptance checks each disk's relative
fixed-point residual separately; a strong rotor cannot mask the weak rotor's
error. The rationalized momentum root retains weak-loading accuracy.
Ordinary and CG force entry points supply a uniform ambient field to this
same local-velocity owner. Invalid or reversed local flow and below-ground surface
samples refuse through the existing envelope receipt without a CG-flow fallback.

Distributed runs select V4 checkpoint integrity and a distinct run identity;
the complete OU state, wind specification and inertial velocity retain their
existing encoding. V3 explicitly selects CG loading and V2 selects legacy
loading. Changing a distributed payload to V3 cannot silently downgrade its
scenario identity. Focused coupling, force and lifecycle regressions cover
independent rotor inflow, the momentum/wash relation, canard-only gust pitch
loads, local-force reconstruction, single pitch-motion ownership and replay.

This remains quasi-steady strip aerodynamics with the existing spanwise prop-wash
overlap approximation, a common straight wake and
reference chord-row circulation ratios, axial disk-center BEMT, reference-flow
parasite/classical induced-drag estimates and reduced longitudinal motion.
It does not supply chordwise gust memory, crossflow dynamics, disk harmonics,
or a validation of historical turbulent flight.

### Recording and scrubbing the selected wind model

Material-bound spatial constructors (`init_with_material_spatial_wind` and
`init_with_material_distributed_spatial_wind`) admit the existing complete
`fs-material::GasState` and replace scalar scenario density **before** trim.
The original trim, wing/canard loads, propellers, rotor dynamics and lateral
loads consequently share its density. `scenario()` exposes this canonical
scenario; `sample_material_air` returns the same accepted wind and gas, with
dynamic-pressure and Reynolds diagnostics relative to an explicit observer
velocity. The complete thermodynamic and transport tuple remains available
without re-evaluation or a second reference-air assumption.

Material runs use V5 checkpoints, retaining the whole numerical gas and an
explicit CG/distributed sampling mode. Every coefficient binds tick-zero and
run identity, including properties not consumed by the current force laws.
The material replay recorders preserve the canonical density and original
initial checkpoint. A contradictory density at restore refuses. Historical
V2/V3/V4 bytes are unchanged. Focused G0/G3/G5 tests in
`simloop/spatial_wind/material/tests.rs` cover actual pressure-dependent flight,
native force/trim equivalence, q/Re, observational noninterference and replay.

The material gas is homogeneous and prescribed: no heat/composition evolution,
altitude stratification or fluid/structure thermodynamics is introduced. The
reference polars remain Reynolds-independent and incompressible; retaining
viscosity, temperature and sound speed does not claim those corrections.
Numerical checkpoints do not serialize material-source evidence or qualify a
gas composition for historical air; that authority stays with the producing
model. The existing flight density envelope still applies.

`replayenv::{record_replay_with_spatial_wind,
record_replay_with_distributed_spatial_wind}` retain the exact V3/V4 tick-zero
checkpoint in the existing envelope layout. Legacy `record_replay` keeps its
checkpoint cadence and serialized layout. `scrub_to_tick(0)` reads the actual
initial state without advancing either physics or the wind; positive ticks
resume the nearest strictly earlier checkpoint and reproduce that step's output.
Both tick-zero digest and run intent must match. Spatial replay requires its
original initial checkpoint even when a later checkpoint is available. Missing,
downgraded, duplicate or wrongly indexed initialization refuses instead of
selecting legacy gusts. The existing deterministic-pilot recording restriction
remains; this does not add an applied human-input trace. The already emitted
damage-model-unavailable terminal code also round-trips through the parser.

## Error model

Typed `Refusal { code, message, ranked_repairs }`. Codes:
`non-finite-input`, `span-outside-domain`, `chord-outside-domain`,
`camber-outside-domain`, `area-inconsistent`, `hinge-axis-outside-prior`,
`pilot-mass-outside-domain`, `component-count-exceeded`,
`component-mass-invalid`, `mass-spec-mismatch`, `hinge-ratio-invalid`.
Caps tested at cap AND cap+1 (workspace law).

The six-DoF spine additionally reports `orientation-outside-group` before
invoking a load callback, `non-finite-loads` for a non-finite force or moment,
`non-finite-state` when finite admitted inputs overflow during a derived step,
and `lie-step-refused` if the canonical integration owner refuses a derived
step. Non-finite start/end time is `non-finite-input` and refuses before load
callback effects.

## Determinism class

Deterministic: pure arithmetic + fs-blake3 digests; the reference design
digest is pinned (golden-bump protocol; one bump recorded when travel_rad
moved from the 0.5236 literal to exact π/6). Spine digests are computed from
canonical output states, so quaternion representatives `q` and `-q` describe
the same callback inputs, trajectory, and digest.

## Cancellation behavior

Synchronous pure functions; nothing to cancel.

## Unsafe boundary

Workspace `deny(unsafe_code)`; no unsafe.

## Feature flags

None.

## Conformance tests

`tests/design_battery.rs`: mass/CG pinned to dossier values (gross 340.2 kg,
CG 29.7% chord); inertia reproduction inside the ±15% band vs the Jex-Culick
lineage with slender-biplane ordering; admission caps at cap AND cap+1
(pilot both edges, fabricated area, hinge prior, hidden-component
mass-spec falsifier, NaN/zero-span); derived panel pinned against
independent hand calculations (canard volume, naive NP, margin arithmetic,
zero gradient at x_h = ¼, self-driving sign, free-worse-than-fixed) plus
the documented naive-vs-Culick NP comparison; digest determinism +
sensitivity + pinned golden. JSONL receipts per case.

## No-claim boundaries

- The published inertias are a single-lineage reconstruction; agreement
  inside the band is a CROSS-CHECK of our decomposition, not a validation
  of either lineage.
- The derived panel's aerodynamic formulas are documented simple models
  (two-surface, no interference, thin plate). Their divergence from the
  Culick neutral point (3.9% c) is recorded data. Real stability claims are
  V-02a territory (E4.6a) and need the full aero model.
- The hinge-moment gradient is a SIGN/SHAPE device over the E1.5 prior; its
  quantitative level inherits the Estimated ceiling (A7a promotion path).
- The lateral state uses `ReducedAeroelasticWarp` and
  `ReducedLateralBuildUp`, with declared inertia and rudder-moment constants.
  It supports bank/heading presentation and reduced control-response claims;
  it is not a six-DOF structural-margin, post-stall, or calibrated-history
  claim.
- No aerodynamic force model lives here; sections are fs-airfoil, planform
  effects are fs-wing.


## Transient wake, shaft and aircraft propulsion

Steady `propcoupling::CoupledStep::power` carries both original BEMT disk
power splits through the accepted coupling evaluation, including spatial
disk adapters and `GearedPropulsion` results. Shaft power agrees with the
physical aerodynamic torque and absolute shaft rate for either blade
handedness. Thrust power uses actual local disk incident air, including the
wing-induced velocity; it must not be mistaken for complete aircraft thrust
work or a closed wing/wake energy balance. Annular induced transfer and
profile-drag dissipation retain their native meanings rather than becoming
additional body heating. The existing geared-propulsion battery checks the
signed shaft-power ports of counterrotating rotors.

`spine::rotors::propulsion::dynamic::DynamicPropulsion` couples two original
uniform axial inflow states to the existing physical geared rotor/aircraft
assembly. `DynamicState` retains both disk-induced velocities alongside native
six-DOF motion and canonical shaft momentum. The apparent-air-mass model is
owned by `fs-airscrew::dynamic_inflow`; its density remains fixed throughout
a trajectory and its air mass is excluded from aircraft inertia.

Each actual aerodynamic evaluation solves the original wing under the
current vector wash, projects incident, wing-induced and other-disk velocities onto
the declared disk axes, and evaluates original blade forces at the stored
induction. The shared carrier-wrench owner applies thrust arms, counterrotating
shaft loads, internal motor torque, and fixed/canted geometry exactly once.
The quasi-steady and dynamic models reuse the same section forces and
body/shaft wrench assembly while retaining distinct public result types.

The initial real loads predict aircraft, shaft and both wake states to
half time. A second real evaluation supplies the accepted native gear map
and the two implicit-midpoint wake-flux steps. The complete result is published
only after every motion and wake stage succeeds; borrowed input state is
unchanged on late refusal. Callback side effects remain caller-owned.
World aircraft/rotor impulse balances and axial reduced-wake impulse/work
balances remain separately inspectable, with no resolved-wake world balance
or whole-fluid/body energy theorem inferred from them.

Four new native coupled cases establish canted startup momentum accounting,
state-dependent motor/gust/body/shaft/wake refinement, exact shaft-sign
convention invariance, and late density refusal. Together with six dynamic
inflow and 37 original BEMT/station/geared/spatial/coupling cases, 47 focused
native cases pass and legacy goldens remain unchanged. Smooth admitted
coupled refinement is second order. The tier retains normal axial,
constant-density, disk-average flow, quasi-steady blade/wing sections and
the existing wash window; it supplies no skew/crossflow loading, swirl,
harmonics, transported wake geometry, blade memory, gimbal motion, historical
calibration, browser wiring or cross-platform qualification.


## Oblique propeller loads in coupled aircraft motion

`spine::rotors::propulsion::dynamic::DynamicPropulsion::new_oblique` explicitly
selects native revolution-averaged oblique blade loads over the existing geared
aircraft and uniform axial wake. The normal constructor keeps the original
axial model. Each actual aerodynamic stage supplies hub incident/wing-induced
velocity, the carrier angular rate and the pointwise other-disk wake to the
native blade owner. Fixed body-frame disk references and physical spin directions are
prepared independently of the signs chosen for mechanical shaft coordinates.

The shared aerodynamic assembly applies the full force at the disk center,
its moment arm about the aircraft COM, and the transverse hub moment. Axial
aerodynamic torque goes only to the original physical shaft port. Consequently
shaft reaction, carrier tilt, transverse loading and aircraft world momentum
use the same real force evaluation without duplicated shaft torque. Both start
and predicted-midpoint stages include oblique loads before advancing the
original aircraft/gear/wake state. No independent mechanics or wing solver is
introduced. The section power account includes hub translation, tilt and
absolute shaft spin; relative shaft power and carrier axial rotation sum to
that same physical absolute-spin work.

The wake remains the named constant-density uniform AXIAL apparent-air-mass
model and its straight axial wash. The scalar momentum/energy receipts are not
a full oblique-fluid impulse or passivity theorem. Ambient/wing incident flow
defaults to hub sampling; the explicit `BladeIncident` policy below can resolve
both fields at native blade sections. Actual carrier tilt and the other disk's
finite wake already vary over the blades. Skew-wake geometry/induction,
blade-phase harmonics, flapping and
unsteady airfoil memory remain outside this model. The existing normal wake
and new section-domain refusals both apply; no negative-speed clipping or
reversed-flow fallback occurs.

Three focused native cases in `tests/oblique_propulsion.rs` reconstruct the
actual carrier wrench difference from blade forces/hub moments, advance a
canted two-propeller aircraft with independent physical world momentum
reconstruction, and reverse all mechanical shaft-axis sign conventions while
retaining an identical physical trajectory. The existing dynamic, geared,
spatial and coupled aerodynamic suites also remain applicable. No historical,
full-workspace, browser, real-time or performance qualification is implied.

## First harmonic disk wakes in native aircraft and shaft motion

`spine::rotors::propulsion::dynamic::first_harmonic::HarmonicPropulsion`
explicitly selects a normal-transport mean-plus-two-mode wake for each
physical propeller. `HarmonicState` retains all six induction velocities
alongside the original carrier pose, twist and generalized shaft momentum.
`disk_frames` exposes the prepared body-fixed mode references so callers
can specify the actual physical field. Mechanical shaft-axis sign choices
do not change these mode directions or their physical forcing.

One original wing solve per actual aerodynamic stage now consumes the
local modal induction at each real wing-strip location inside the existing
straight axial wash window. The same field reaches every native blade
quadrature point. Blade thrust and actual transverse hub moments drive the
three modal air-inertia ports; all body forces, lever arms, transverse hub
moments and aerodynamic shaft torques pass through the original shared
wrench and geared mechanics owners exactly once. No duplicate section,
wing, mechanical or nonlinear wake solver is introduced.

The initial force stage predicts aircraft, shaft and all wake modes to
half time; the actual predicted-midpoint stage advances the original
mechanics and each implicit transport step. The wake input work and the
section induced work use their actual respective midpoint states.
`coupling_work_residual_j` explicitly exposes their difference rather than
calling a partitioned quadrature defect dissipation or silently claiming
coupled passivity. The modal frame-rate term uses carrier rotation about
the disk normal, never physical shaft spin. Original immutable-state and
late-failure atomicity remain; callback effects are caller-owned.

The selected wake closure is still a normal-flow, zero-skew transport
approximation even when the blade law resolves incident crossflow. The
existing straight wash geometry is a finite-window transfer approximation;
it is not a free wake, skewed induction solution, full world-fluid momentum
account, wake-distortion model or total fluid/body energy theorem. Incident
air/wing-induced flow defaults to each disk center, with the explicit
`BladeIncident` policy below resolving spatial variation at actual blade
sections. Carrier tilt and other-disk wake loading retain their original
pointwise owners. Finite blade phase, higher modes, flapping, cyclic
pitch, swirl and wake deformation remain unimplemented in this tier.
Apparent air inertia remains excluded from
the aircraft's rigid inertia.

Four native cases in `tests/harmonic_propulsion.rs` establish zero-mode
parity with original oblique aircraft loads; nonuniform induction reaching
real blades and wing strips; independent physical aircraft/shaft world
momentum reconstruction over 24 canted-rotor ticks; full-state and coupling
work refinement; mechanical-axis sign parity; and late-density refusal with
deterministic retry. Measured 24-tick residuals are `8.73145e-13 N s` linear,
`1.58448e-11 N m s` angular, `2.94253e-15 N s` modal and `2.10943e-15 J`
within-wake storage/work. At 24/48/96 steps, the complete-state error against
384 steps is `3.56045e-4`, `8.81038e-5`, `2.09938e-5`; accumulated absolute
induced coupling-work defects are `1.45322e-3`, `3.65155e-4`, `9.15029e-5 J`.
The state error is a numerical fixture norm, not a common physical unit or
trajectory certificate. All ten new and 55 original propulsion/airscrew
cases pass in the focused actual-source native build. This is not a full
Cargo/DSR, cross-platform, browser, real-time or physical validation gate.

## Skew-flow harmonic wakes in actual turning aircraft motion

`HarmonicPropulsion::new_skewed` selects the full Pitt--Peters first harmonic
transport from the `fs-airscrew` source owner. `new` retains its original
normal-transport choice. The dimensional source gain, moment signs and
normal-throughflow limits are specified in the airscrew contract and Chen,
NASA TM-102219, equations (13b), (14), (17), (18):
<https://ntrs.nasa.gov/api/citations/19900006622/downloads/19900006622.pdf>.
This extends the earlier explicitly zero-skew transport assumption; it does
not replace the shared section, wing or geared mechanical owners.

Each actual aerodynamic stage projects the same wing/ambient incident
velocity into its fixed physical disk frame. The current transverse
direction changes only the source operator; stored wake modes retain their
body-fixed meaning. `disk_incident_mps` and optional `skew_transport` expose
these actual stage inputs and coefficients. Mean and both harmonics reach
the original radial/azimuth blade points and existing wing-strip wash
window. Their actual thrust and transverse hub moments drive the modal
ports. Rotor body/shaft forces, torque, carrier reaction, lever arms and
motor work still use the original mechanical update exactly once.

The optional `skew_inflow` entries contain the same accepted modal steps
also exposed in `inflow`; they are additional accounting detail, not extra
updates. Source cross work is signed and does not become invented heat.
The original `coupling_work_residual_j` still measures implicit wake input
work minus actual predicted-midpoint blade induced work. A late wake,
blade, wing or input refusal returns no changed borrowed state; callbacks
remain caller-owned. Source iteration caps and tolerances are explicit.

Three additional actual-source aircraft cases exercise a time-varying
crosswind and turning disk velocity; exact shaft-axis/ratio sign parity;
unchanged aerodynamic loads at identical induction but changed derivatives;
independent physical world momentum reconstruction; signed source work;
full-state and partitioned-work refinement; exhausted iteration budgets;
and late refusal/retry. A 24-tick, 96 ms canted-rotor trajectory measures
`2.67193e-13 N s` linear and `1.84520e-11 N m s` angular impulse defects,
with largest modal defect `4.76659e-11 N s`, within-wake energy defect
`1.89884e-11 J`, and at most four coefficient iterations. Its accumulated
signed cross transport work is `-2.21195e-3 J`.

The complete coupled refinement check uses the initial smooth 6 ms interval
of that fixture: 8/16/32 steps versus a 256-step reference give state-norm
errors `1.79455e-4`, `4.49075e-5`, `1.10028e-5` and accumulated absolute
coupling-work defects `5.12082e-5`, `1.22387e-5`, `2.97900e-6 J`. This norm
mixes native coordinates and is a numerical convergence measure, not a
common physical unit. Longer refined versions encounter nonuniform
convergence or the existing quasi-steady wing solver's 600-iteration
refusal; no wing tolerance, cap or closure was relaxed. The 96 ms momentum
ledger is therefore not a convergence certificate for that entire flight.
Independent source-level nonlinear refinement remains separately tested.

All nine new and 65 prior airscrew/propulsion cases pass in the focused
pinned native build (74 cases across twelve suites). This does not claim
full Cargo/DSR, cross-platform, real-time or physical validation. Existing
positive section throughflow/angle domains, disk-center incident sampling,
straight finite-window wash geometry and constant air density still apply.
The extension supplies the skew transport matrix; it does not supply a
transported/skewed far-wake geometry, resolved intra-disk wake loading,
intradisk ambient wind gradients, reverse-flow regimes, higher modes,
swirl, flapping/cyclic pitch, or a complete fluid/body energy/momentum model.

## Direct momentum transfer between propeller disks

Steady, uniform-dynamic and first-harmonic propulsion now include the other
propeller's current wash in the incident velocity used by the original blade
forces and wake momentum law. Self-wash is excluded from this transfer: each
disk's native blade/wake law still owns its own induction exactly once. The
wing-induced field is added once from the same accepted wing circulation.
Steady transfer uses the current joint coupling iterate; transient transfer
uses both physical induction states at the actual force stage.

The receiving-disk mean integrates the existing straight finite-2R-window
field over its physical area. Parallel equal-radius disks use exact circular
overlap and its exact centroid, so partial overlap remains live even when
neither center is inside the other wake. Canted disks use a deterministic
32-annulus by 64-azimuth equal-area midpoint rule. The prepared coefficients
include transverse first moments; source harmonic modes are evaluated in
their own physical disk frame before their receiving-area mean is transferred.
Geometry is prepared outside steady iteration and once at dynamic/harmonic
model construction; no additional nonlinear wing solve is introduced.

This extends the declared actuator field, not a resolved mutual wake. Axial-only
blade models consume one mean external wash vector. Oblique and harmonic
blades additionally resolve the same finite field at every section as described
below. Wake deformation/contraction, swirl and mutual-wake world energy
closure remain outside this tier. The canted area rule has a focused
refinement test, not a uniform quadrature-error certificate. The existing
wing-strip wash model and gates remain unchanged. Tests cover analytic lens
area/centroid, canted refinement, separated/coaxial/partial-overlap disk
inflow, native blade/momentum consumers, and zero-harmonic parity.

### Pointwise external wake loads on oblique blades

`DynamicPropulsion::new_oblique` and all `HarmonicPropulsion` variants use the
other rotor's physical mean/first-harmonic wake at every original radial and
azimuth blade sample. Partial overlap therefore changes local force and creates
the corresponding asymmetric hub moment, even when neither hub lies inside
the source cylinder. Each source's harmonic coordinates use its own fixed
disk frame. Self-wash stays in the original blade owner; the receiving-area
mean is subtracted from the prepared variation before adding it to the same
mean incident velocity used by transport. Neither contribution is doubled.

Fixed geometry, cylindrical support, axial influence and source-mode
projections are prepared once at construction in two bounded arrays. The force
loop needs only three coefficient products per section and no geometry queries,
trigonometry or allocations beyond the original blade quadrature. Uniform
axial-only models do not allocate those arrays. Zero external modes retain the
original native blade path. The accepted force, transverse hub moment and
separate signed shaft torque pass through the original shared aircraft wrench
and geared midpoint mechanics, including continued and skewed harmonic stages.

The native blade result reports `incident_variation_work_w = sum(dF dot delta_V)`
separately from mean-relative hub work, own induced work and profile heating.
That signed spatial-air term closes the local section power identity; it is not
extra carrier work or a total two-wake energy balance. The scalar/modal wake
transport still consumes its declared receiving-disk mean, not a newly resolved
transport PDE. Ambient/wing-field gradients use the separate opt-in policy
below. No wake deformation, contraction, swirl or full fluid/body conservation
theorem is introduced.

The partial-overlap harmonic test independently reconstructs pointwise support
and source coordinates, checks the new hub moment and work, and retains exact
zero-harmonic parity and absence of self-wash. Canted force/arm/shaft and actual
trajectory momentum regressions exercise overlapping disks. Native blade tests
independently integrate the original section law and spatial-air work.

### Spatial wind and accepted wing fields at real blade points

`dynamic::BladeIncident` is an explicit force-resolution choice for
`DynamicPropulsion::new_oblique` and every `HarmonicPropulsion` variant.
`with_blade_incident` selects `Hub` (the unchanged default), `WindField`, or
`WindAndWingField`. Selecting either field mode on the axial-only dynamic
constructor refuses with `dynamic-propulsion-invalid`; no hidden quadrature or
different force owner is substituted. The same selected path is used at start,
midpoint, skewed and continued stages.

The wind callback keeps its original contract: it returns world-frame material
air velocity [m/s] at an actual world position [m]. For a body point `p`, the
stage samples `x_world + R_body_to_world * (p - cg_body)`, at that stage's time
and pose. It caches the admitted pose and each already-required hub wind
sample. `WindField` then supplies only
`R_world_to_body * (wind_world(hub) - wind_world(blade))` as the additional
relative-air velocity. The subtraction occurs before rotation and before any
large carrier velocity is introduced. Native carrier tilt, absolute shaft
rotation and own induction remain owned by the original blade evaluator;
none is resampled or doubled. Spatially uniform world wind yields exact zero
increments and retains the original loads, including partial mutual overlap.

`WindAndWingField` additionally evaluates the original free-air wing kernel at
each blade point using the same accepted circulation and stage freestream.
Subtracting its already-included hub value adds exactly the local wing-field
variation. The original kernel/core treatment and one wing solve per force
stage are retained. No nonlinear wing solve runs inside the blade loop. Neither
mode allocates per blade sample. Wind mode adds one material-air callback per
sample; wing mode also adds an O(panel-count) native kernel evaluation. The
already-prepared mutual-wake geometry and section quadrature are reused.

All added fields enter the original section force law. The resulting full
force, transverse hub moment and separate signed shaft torque reach the same
aircraft/gear owners. Their combined signed material-air work remains in
`incident_variation_work_w`; it is neither carrier work nor dissipation. The
scalar/modal wake law continues to use the original hub ambient/wing velocity
and receiving-disk mean mutual wake. Resolving blade loading does not claim a
resolved wake transport PDE, full fluid momentum account, blade-phase model
or closed wing/wake energy balance.

World-position, material-velocity and transformed-increment failures refuse
before publishing a load. Callback refusals keep their original code, message
and repairs through the native blade adapter. A refusal at a later blade or
predicted-midpoint stage leaves caller-owned physical state and continuation
anchors unchanged; callback side effects remain caller-owned.

The incident-field cases in `tests/harmonic_propulsion/incident_field.rs` check
exact uniform-wind and zero-mode parity, independently reconstruct native
blade loads and signed spatial-air work under world shear and wing variation,
and verify actual carrier wrench assembly. They exercise both canted physical
spin directions, equivalent mechanical sign conventions, rotated world poses,
late blade/midpoint refusal with exact retry, and actual world momentum along
a trajectory with evolving harmonic wakes.

## Actual stage-by-stage wing continuation in skewed flight

`first_harmonic::ContinuedHarmonicState` pairs the existing complete physical
`HarmonicState` with an explicit `fs-wing::nonlinear::continuation::WingContinuation`.
`initialize_continuation`, `evaluate_continued` and `step_continued` select the new
wing continuation lane. The original entry points retain their original cold
wing solve. The actual initial force stage continues from the caller's anchor;
its accepted circulation seeds the actual midpoint force stage. The returned
physical endpoint carries the accepted midpoint circulation as the next initial
guess. Both stages resolve current section equations and apply current loads;
neither force is lagged. The original aircraft/shaft mechanics, blade section
quadrature, skew-inflow owner, force/moment arms and impulse/work accounts remain
the only physical owners. The anchor adds no aerodynamic momentum state.

The exposed long-flight failure had two distinct causes. On the original
`vel_z=0.8 m/s` positive-gust trajectory, the old Picard iteration refuses at
8 ms after its 600-attempt cap even while the section solution is attached and
the row split is well conditioned. Caller-owned warm starts alone did not cure
it. The new native Newton--Krylov continuation crosses that failure at the same
`1e-10` actual componentwise section gate. Continuing the same positive-gust
trajectory later approaches cancellation of a signed linear chord-row total;
the declared row-condition bound correctly refuses. The original
`vel_z=1.5 m/s` fixture also has near-filament and very large effective-incidence
roots; the admitted continuation API explicitly refuses it. A low algebraic
residual on that legacy fixture is not evidence of realistic long-trajectory
accuracy.

The full coupled refinement regression uses the original attached
`vel_z=0.8 m/s` initial state and the opposite smooth crosswind turn,
`wind(t,p)=[0.02*p_y, -0.5*sin(20*t), 0.1*cos(20*t)] m/s`. It completes the actual
96 ms interval at 24, 48, 96 and 384 physical steps with canted counterrotating
disks, moving attitude, local wing velocities, first-harmonic skew transport and
state-dependent aircraft/shaft loads. The section equations at every real
initial and midpoint stage are independently reconstructed from the public full
horseshoe field minus an independent finite-bound field. Their maximum measured
normalized residual is `3.069e-11`; all effective incidences remain attached.
The run does not weaken the force residual or shorten the physical interval.

Against the 384-step reference, the existing full-state distance errors at
24/48/96 steps are `3.4904523143e-4`, `8.6389460900e-5`, `2.0587886173e-5`.
The accumulated absolute induced-work coupling defects are `1.4235880765e-3`,
`3.5729731447e-4`, `8.9494495170e-5 J`. Both refine at the expected second order.
Across those runs and the reference, independent solid-body/rotor world linear
momentum defects stay below `7.48e-12 N s`, angular momentum defects below
`9.03e-11 N m s`, reduced-wake impulse defects below `6.92e-11 N s`, and wake
energy-balance defects below `3.59e-11 J`. No stage uses more than six Newton
attempts or 44 Jacobian products. Reversing equivalent shaft-axis/ratio/blade
handedness conventions gives exactly the same complete continued run and work
account. Late source and exhausted-Newton refusals leave both physical state and
wing anchor unchanged; exact replay is tested.

These measured bounds concern this admitted 96 ms convergence fixture. They do
not certify arbitrary maneuver duration, separated-wing roots, a moving far
wake, unsteady section hysteresis or full fluid/body momentum closure. The
previous explicit induced-work partition defect remains a quadrature defect,
not dissipation or a passivity guarantee. Row-cancellation, incidence and
filament domains remain explicit physical/numerical admission limits rather
than hidden regularization.

The focused final sweep passes 107 cases: 30 wing/continuation cases and 77
airscrew/propulsion cases, including the ten-case harmonic flight battery.
The native harness compiles current wing and propulsion sources and the actual
Newton/FGMRES/Armijo numerics; it omits the solver's two `fs_exec::Cx` convenience
wrappers and unrelated aircraft surfaces. This does not substitute for a full
workspace Cargo/clippy/CI gate, which was not run in this environment.

## Gravity in longitudinal perturbations

`longitudinal::linearize` holds pitch fixed while perturbing body air
velocities. The original `aircraft::force_buildup` retains its level-flight
trim convention, where pitch equals incidence. Its inserted weight is removed
before taking velocity derivatives; actual attitude-dependent gravity enters
only the pitch column. Angle-of-attack changes therefore cannot create
spurious gravity terms in the speed or incidence derivatives. Aircraft force
assembly, linearization, airborne motion and rail force projection share the
same terrestrial gravity constant and weight-removal operation.

The focused G1/G3 checks distinguish incidence from attitude using an exact
gravity-only manufactured field and compare full native aerodynamic force
derivatives after changing mass at a fixed operating point. The latter must
remain mass-independent before division by inertial mass. This corrects the
linearized poles and all augmented/pilot models consuming that matrix; it does
not change the underlying aerodynamic solve or the original level-flight trim.
The historical model still uses uniform standard Earth gravity. Spatial
planetary gravity is not inferred from this correction.

The focused x86 native sweep passes 27 cases: 11 motion tests, seven V-02a
tests, seven augmented-spectrum tests, and the actual assistance/pilot flight
behavior regressions. The original longitudinal source reproduces the previous
V-02a digest and fails the new mass-scaling regression. Correcting its gravity
derivatives changes nominal doubling time from 0.6204376902 to 0.6670863177 s;
all physical acceptance thresholds remain unchanged. Only measured x86 golden
pins were updated; other ISA pins still require native remeasurement.

This source-selected build uses the original aircraft/wing/propeller and
Newton/FGMRES numerics, excluding only unrelated module surfaces and the two
unused solver `Cx` convenience wrappers. It does not establish a full Cargo or
DSR gate; full-runtime builds were interrupted while compiling asupersync.


## Prescribed material atmosphere in the evolving flight lifecycle

`SimLoop::init_with_profile_spatial_wind` and its distributed counterpart bind
an immutable `MaterialAirProfile` to the original native trajectory. Initial
trim and prelaunch state use the source gas at the 2.4 m rail height. Each
beginning/midpoint load stage samples its own CG height, and the same density
drives native wing/canard forces, propeller thrust and torque, aeroelastic
dynamic pressure, and reduced/distributed lateral response. The value in
`scenario().rho_kg_m3` remains the initial rail density. `sample_material_air`
returns the complete original gas at the requested height without advancing
wind or trajectory. Predicted endpoints are profile-admitted before publishing
motion and its corresponding warm start.

V6 checkpoints retain the profile definition, geometric ground datum, sampling
mode and complete existing wind/flight state. Uniform profiles retain all twelve
gas fields; the standard atmosphere retains its model discriminator and ground
altitude. Profile identity binds into tick-zero and run identities. Restore
admits the initial density and saved height before publishing state. Matching
replay recorders retain the initial V6 checkpoint; existing V2–V5 modes keep
their layouts and homogeneous paths retain their existing arithmetic.

The original material owner supplies hydrostatic pressure and coherent gas
thermodynamics/transport. Density is frozen across the aircraft within each
load stage; distributed mode resolves incident velocities across strips and
disks. Existing reduced flight dynamics, constant flight gravity, engine torque
law and Reynolds-independent polar remain in force. This adds no weather
evolution, compressible/anelastic atmospheric dynamics, buoyancy, pressure work
or high-altitude aircraft validation claim. A prescribed solenoidal velocity
field with stratified density does not imply divergence-free mass flux.

Six focused regressions cover native trajectory/rotor changes with altitude,
stage-density coherence through midpoint/warp/lateral loads, homogeneous
equivalence, profile-domain refusal before motion publication, complete source
encoding, and checkpoint/replay continuation. Source/parser review is distinct
from native execution, which remains unverified for this addition.

## Coupled rotor and aircraft momentum in the augmented model

`augmented::build_engine` takes rotor torque derivatives from the accepted
`OpenLoopDesign::force_buildup` loads in body axial/transverse velocity, pitch
rate, canard deflection and common propeller speed. For the existing equal
per-prop inertia `I_r`, the shared-speed equation is
`2 I_r omega_dot = 2 Q_engine - Q_left - Q_right`; differentiating the
aircraft's mean torque imbalance and dividing by `I_r` preserves that total
shaft angular-momentum balance. Both hub velocities and the converged wing
induction therefore affect rotor acceleration.

The same speed-perturbation load pair supplies the complete rotor-to-body
column: axial and vertical forces and pitch moment, including prop wash and
each propeller's actual moment arm. The supplied rigid matrix and canard
control column retain their owners unchanged. Attitude is independent of air
velocity; gravity cancels from the speed column at fixed incidence. The model
still has one synchronized rotor-speed state, quasi-steady aerodynamic loads
and a longitudinal rigid block. It does not add independent rotor modes or
six-DOF gyroscopic dynamics.

The native numeric-source slice passes 28 focused tests: 12 aircraft, seven
longitudinal and nine augmented cases. Independent off-trim torque/force
samples use five-times smaller perturbations than the engine; unequal hub
arms and disk relabeling expose missing load paths. Both new coupling tests
fail against the previous augmented owner, whose existing x86_64 golden is
reproduced before the correction. Only the measured augmented x86_64 pin
changes; aircraft and longitudinal pins and all physical bands are unchanged.
As in the preceding velocity checks, this bounded source-selected build
omits the unused continuation branch and unrelated runtime modules; it is not
a full workspace Cargo or lifecycle validation.

## Steady angular wake momentum in two-way geared propulsion

`propcoupling::swirl::coupled_step` selects the native two-component
`fs_airscrew::swirl` solver inside the original two-way wing/rotor iteration.
Actual swirl changes the blade thrust, required shaft torque, wing wash and
accepted circulation. `CoupledSwirlStep` retains both complete annular
momentum solutions from the same accepted evaluation; no diagnostic rotor
re-solve or one-way fallback supplies the published loads. `CoupledStep<P =
BemtPower>` preserves the existing default while this explicit variant carries
the native axial, angular-wake and profile power components separately.

The original iteration, correction budget, rollback and per-disk residual
gate are shared without arithmetic changes. Rotor continuation uses fixed
caller-supplied annular seeds throughout each outer solve, so rejected outer
trials never become new branch anchors. The complete accepted solution can
produce the next `WarmStart`. These seeds are numerical data, not added air
inertia or physical momentum. Model identity includes the selected swirl
law and its native bounded-solve controls; the ordinary coupling digest stays
unchanged.

`spine::rotors::propulsion::swirl::SwirlPropulsion` binds this selection to an
existing admitted `GearedPropulsion`. It reuses the same material-air sampling,
spatial disk axes, wing force/moment assembly and absolute physical shaft
rates. Required aerodynamic torque enters the existing external shaft ports
with the actual blade handedness. It is not added again to the carrier-only
wrench. Actual start and midpoint wing/swirl evaluations drive the original
geared six-DOF step, including rotor acceleration, internal motor recoil,
gyroscopic response, spin-plus-orbital momentum and physical work.

Five regressions in `tests/swirl_propulsion.rs` establish two-way force and
circulation changes, independent disk-inflow/native-solution reconstruction,
signed shaft power, co-rotating torque ownership, actual global momentum
change versus external impulse, regenerative windmilling acceleration, canted
shaft world-frame covariance, and bounded/late-stage atomic refusal.
Existing no-swirl force/torque and coupling golden arithmetic is preserved.

The spatial wash supplied to wings and other disks remains the existing
axial actuator-disk field. This selection adds no transported helical-wake
velocity, dynamic swirl memory, radial pressure/contraction, reversed-flow
closure, combination with harmonic wake models, historical calibration or
whole-aircraft/fluid conservation theorem. Native original-module checks
are distinct from full-workspace Cargo, cross-platform or real-time gates.


## Impulses in prepared passive shaft motion

Prepared `ShaftMotion` also accepts instantaneous carrier and independent-shaft
impulses through `apply_impulse`. It updates the retained native canonical
momentum and actual Flyer velocity while preserving pose, relative phase, spring
potential and the prepared spring/drag model for subsequent continuous steps.
The event retains native actual-impulse and kinetic-increment receipts and
separate motor, aerodynamic and carrier work. Original Flyer transfer checks
and late-refusal atomicity apply. No collision/restitution law, shaft-lock mode,
time advance or event-to-heat allocation is inferred.

## Opt-in torque transport into downstream wing and carrier momentum

`propcoupling::swirl::wake::coupled_step` projects the actual accepted native
annular blade torque into a downstream azimuthal velocity before solving wing
circulation and assembling its forces. Each original radial cell carries
`L_i = Q_i/(rho U A_i)` and `v_theta = L_i/r`, so independently integrating
`rho U r v_theta` over its physical area recovers the original signed `Q_i`.
Prandtl factors are not inserted twice. Native at-disk induction `u` is distinct
from downstream swirl; only a no-loss, matching-flow narrow-annulus limit
recovers `2u`. Physical spin relative to the positive-thrust axis determines
the direction, and material-relative wing velocity receives the opposite of
the physical fluid swirl. Signed windmilling reverses angular transport.

A bounded relaxed field/wing/rotor fixed point surrounds the unchanged inner
two-way axial coupling. The same prepared geometry/operator and actual native
disk-inflow samples are reused. Caller rotor branch seeds remain fixed through
all outer trials. `CoupledWakeStep` retains the exact swirl velocity used in
the accepted wing force solve, the native source fields from that evaluation,
and the maximum used-versus-predicted strip-velocity discrepancy. Exhaustion
refuses; there is no one-way fallback. The native result and its usual inner
coupling residual remain available separately.

The opt-in `spine::rotors::propulsion::swirl::wake::WakePropulsion` evaluates
these loads at both initial and midpoint mechanical stages. It supports the
admitted local section frames and derives spin from mounted physical shaft
axis plus blade handedness. Actual wing forces/moments and signed native
shaft torque reach the existing geared carrier/rotor integrator once each.
The ordinary `SwirlPropulsion`, axial coupling, historical presets and their
default payload types retain their existing behavior.

The projection deliberately retains the original fixed-cylinder axial wash
approximation. Its changing downstream axial speed changes mass flux at
fixed area: no contraction/entrainment solution, radial pressure balance,
free/helical wake or whole fluid/body conservation theorem is supplied.
These are unobstructed rotor source fields, not a reconstructed exhaust after
a wing has extracted momentum. The observation window is strictly
`0 < downstream distance < 2R`. In-annulus queries on the disk or beyond that
window refuse; upstream, lateral and hub samples are zero. Angular momentum
is not silently discarded at a numerical cutoff.

Resolved swirl kinetic power is integrated analytically, including its finite
radial-cell logarithm, and compared with each native `(u/r)Q`. The sum of
absolute cell discrepancies is bounded by caller-admitted absolute watts plus
a relative allowance against native total swirl power. Monotone axial speed
and convexity in `1/U` make the two inside-window endpoint limits a bound on
the complete window. Opposite cell errors cannot cancel. The signed defect
and changing mass flux remain visible; unresolved angular structure and model
error are never relabeled as heat. The default relative projection allowance
is 0.5 of the swirl channel, with zero absolute allowance; callers can require
a tighter budget. This tolerance is explicitly separate from numerical
convergence. Receiving disks potentially intersecting another downstream or
coplanar source cylinder refuse. A wing sample receiving two nonzero
overlapping angular fields also refuses without a wake-mixing model.

Six end-to-end tests and two analytic projection tests cover actual
handedness-dependent wing force, independent per-annulus polar integration of
torque and resolved kinetic energy, zero/low-swirl behavior and second-order
radial energy-error refinement, shaft-axis conventions, regenerative
windmilling, actual midpoint global linear/angular momentum, world covariance,
source-window limits, power-budget/iteration refusals, unsupported disk/wake
overlap, and late-stage atomicity. In the symmetric downstream-wing fixture,
reversing physical handedness changes vertical aerodynamic force from
`+1.50626572984 N` to `-1.50626572984 N`; axial-only wash gives zero. The
unmodified existing swirl/geared/spatial/coupling batteries pass another
38 cases. These 46 native checks use full selected current source modules,
fresh native airscrew and Lie/gyrostat builds, and original numerical kernels;
the cached solver omits only unused runtime `Cx` wrappers. They are not a
full workspace Cargo, cross-platform or real-time gate.

## Material annular wake transport into flight dynamics

The opt-in `WakeTransport::MaterialStreamtubes`, selected through
`WakePropulsion::with_transport` or `coupled_step_with_transport`, replaces
downstream axial and angular wash with the actual native station distribution.
For each original disk annulus, `m_i = rho*(V+a_i)*A_i` and
`ell_i = Q_i/m_i` remain constant. At downstream distance s the prescribed
axial development is `U_i = V + (1+s/sqrt(s*s+R*R))*a_i`, and its area is
`m_i/(rho*U_i)`. Adjacent areas build contiguous material annuli around an
explicitly unloaded fixed-radius hub. Powered annuli contract; braking annuli
can expand. Original signed blade torque and Prandtl losses are never rescaled.

Changing area also requires radial fluid motion. For inner squared radius b,
the same geometry gives `v_r = (U*b' - U'*(r*r-b))/(2*r)`. This satisfies
cylindrical continuity inside every source annulus and
`v_r/U = dr_boundary/ds` on both material interfaces. Discontinuous axial
velocities may therefore represent tangential sheets without exchanging mass
between adjacent annuli. The material-relative aerodynamic increment is
positive along the thrust axis for positive axial induction and opposite the
physical radial and azimuthal velocities.

The original inner coupling still owns native wing/rotor solves and its branch
seeds. A small private callback retains the actually evaluated axial wash x,
distinct from its next predictor G(x). The outer field correction removes the
actual downstream legacy wash returned with the accepted circulation, retains
original upstream wash, and inserts the complete material-streamtube field
once. The accepted force velocities remain in the original coupled result.
The published `used_strip_swirl_velocity_mps` and warm-start field names are
retained; in material mode they hold the all-component numerical correction.
`strip_velocity_correction_mps()` exposes that meaning without another mutable
velocity owner. Both initial and midpoint geared flight stages use this path.

`material_annuli(s)` returns the physical radii, interface slopes, mass and
angular fluxes, resolved swirl power, radial kinetic power and signed axial
kinetic transport above the native incident speed. Radial kinetic power uses
a stable positive completed-square analytic annular integral; narrow cells
avoid subtracting nearly equal logarithms. `cross_section(s)` reports constant
total material mass flux and native signed angular momentum flux, with its
axial speed interpreted as an area-weighted mean. The fixed-cylinder default
keeps its original numerical behavior and uniform-speed interpretation.

Power admission bounds every cell's resolved swirl discrepancy over the
complete 0..2R interval. Independent positive area/velocity bounds remain
conservative when annuli have mixed signed induction and the cumulative
radius is nonmonotone; opposite errors cannot cancel. The same conservative
outer-radius bound rejects receiver disks reached by an expanding wake.
Unsupported source overlap, reversed throughflow, unresolved annular width,
window escape and nonconvergence still refuse.

This is a steady source-streamtube closure, not a radial-pressure/Euler solve,
a solution for the surrounding bypass flow, a turbulent mixing/entrainment
law, or an exhaust reconstructed after a wing extracts momentum. The model
does not claim complete fluid/body energy or linear-momentum closure.
Resolved meridional kinetic transport is exposed for that future pressure/
energy closure and is never credited as heat or removed from rotor torque.
The original finite observation window and disjoint source-disk restrictions
remain explicit.

## Annular source fields on independent shafts

`spine::rotors::propulsion::passive::wake::WakeShaftPropulsion` connects the
existing native annular angular-momentum field to independent physical shafts.
It shares the prepared `WakeModel` with geared `WakePropulsion`: mounted shaft
axes and blade handedness determine physical swirl signs once, and the same
native wing/rotor/field solve consumes each actual signed axial-momentum/J speed.
Fixed-cylinder transport remains the default; material streamtubes explicitly
include downstream axial replacement, radial motion and contraction/expansion.
The existing field convergence, source separation and power budgets are unchanged.

Both actual mechanical stages use the retained `ShaftMotion` through the common
independent-shaft load callback. Native signed blade torque enters its external
shaft port exactly once; downstream wing force/moment enters the carrier wrench.
Numerical warm starts remain local through the step. A field, wind, source-budget
or late mechanical refusal preserves both canonical mechanical states and caller
seeds. The native owner retains motor work, bearing heat, spring storage and
carrier recoil; this source-field selection adds no mechanical momentum owner.

Focused tests compare exact geared/independent snapshots for both transport
choices and opposite physical shaft-axis conventions, resolve unequal speeds
after impulses with independent spin-plus-orbital momentum and mechanical work
accounts, check native torque and material mass flux downstream, and exercise
signed windmilling plus budget, midpoint and late native refusal/retry. These are
selected-source native checks, not full-workspace, cross-platform or real-time
gates. This integration makes no additional timestep-convergence claim.

The annular source field remains steady and distinct from physical harmonic
wake memory. Optional native wing-branch continuation is described below; it
does not add physical air inertia. The finite observation window, disjoint
source restrictions and explicit projection discrepancy remain in force;
there is no new unsteady wake law or complete fluid/body conservation claim.

## Wing continuation through annular source-field coupling

`propcoupling::swirl::wake::coupled_step_continued` reuses the original
`fs-wing` Newton/FGMRES/Armijo continuation owner inside the existing two-way
rotor and outer annular-field iterations. Every inner and outer trial holds
the caller's complete accepted wing anchor fixed. Numerical wash/rotor seeds
remain separate; a rejected trial never advances the caller's branch.
The private inner result carries its native continuation alongside the actual
nonlinear report, and captures it only at the accepted coupling residual.
Only a complete outer field acceptance returns `ContinuedWakeStep`, pairing
the exact circulation/incidence anchor with its actual used force velocities.
The existing `CoupledWakeStep` type and ordinary entry points remain intact.
Continued coupling identity includes every declared wing-domain/solver limit.

The shared prepared `WakeModel` exposes `evaluate_continued` and
`step_continued` on both geared `WakePropulsion` and independent
`WakeShaftPropulsion`. Evaluation with no previous anchor establishes the
initial complete branch without advancing mechanics. During a step, the
complete accepted initial aerodynamic solution supplies the local midpoint
anchor. Both physical mechanical stages retain their original owners and
signed torque/work semantics. Independent shafts publish the caller's wing
anchor only after native mechanics succeeds; the geared result returns its
midpoint anchor with the accepted state. Native source, field, wing-domain
and late mechanical refusals publish no partial branch or mechanical state.

Four focused regression groups exercise both field choices: changing-load
80 ms nonlinear flight with unequal shafts and actual Newton work; exact
stage circulation/incidence/force-anchor correspondence and independent
accumulated momentum/work balances; fixed caller-branch limits and a case
where continuation succeeds while original cold Newton refuses; atomic
midpoint/late-native refusal and deterministic retry; and geared-stage
replay/parity. The existing 13 annular consumer checks also pass. This is
selected-source native evidence, not a workspace, cross-platform, real-time
or new timestep-order gate. Continuation supplies an algebraic branch anchor,
not dynamic wing stall, wake inertia or an additional fluid closure.

## Prepared tree evaluation of connected wake fields

`PreparedWakeField::prepare_tree` prepares the existing VPM source tree once
under separate byte and work limits. Its returned `PreparedWakeTreeField`
borrows the immutable source identity and reuses source geometry, moments,
source indices and traversal scratch. Physical and reflected ground probes use
one tree. The extra tree byte limit excludes the already owned exact snapshot.

`PreparedWakeTreeField::sample` runs the original field-service assembly,
tensor invariants, masks and finite-difference divergence through this prepared
source evaluator. One request-wide work limit covers all physical/image center,
FD and core queries. A full Jacobian query accepts an expansion only under both
velocity and Frobenius-gradient error limits; exact leaves retain the original
finite-core kernel and zero-core singular masking. Core masks use exact source
distances with conservative tree pruning, so their final pass does not restore
a full source scan at every distant grid point.

The requested absolute limits are shared between physical and image systems.
The `BoundedWakeFieldSample` wrapper retains the requested accuracy and actual
work and maximum error allowances beside the original arrays. Curl and analytic
divergence bounds follow from the Jacobian norm. The six FD velocity probes use
a tightened tolerance accounting for division by the original FD step; their
reported propagated approximation allowance fits the requested gradient limit.
FD discretization bias, atmospheric model error, nonlinear Q/lambda2 uncertainty,
trajectory error and directed-rounding certification remain outside that claim.
Invalid source points retain the existing invalid mask. Zero tolerance evaluates
original source order exactly, with compensated sums; floating-point reduction
may differ from the older un-compensated exact field by roundoff.

Focused regressions compare exact-mode fields and source identity, distant
accepted expansions and achieved error bounds, image axial-curl reflection,
exact core masks, request-wide work exhaustion, cancellation/retry and zero-core
masking. Native Rust compilation and execution were unavailable for this
extension; the independent numerical expansion cross-check is separate evidence.
