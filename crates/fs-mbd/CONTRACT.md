# CONTRACT: fs-mbd

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
- `Wrench` carries a constant world-frame force and principal-body-frame torque
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

## Invariants

- Public construction rejects non-finite mass-property, state, gravity, wrench,
  duration, and quaternion inputs at the relevant boundary. State, pose, and
  gravity fields cannot be forged externally; public integrator operations also
  revalidate their inputs before producing a receipt or diagnostic.
- The public orientation remains normalized and uses a canonical quaternion
  sign after every successful attitude update.
- The declared force and gravity are world-frame quantities; the declared
  torque and angular momentum are principal-body-frame quantities. Diagnostics
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
- No holonomic loop constraints, nonholonomic rolling constraints, RATTLE
  projection, generalized-alpha lane, general Newton/Krylov API, or constraint
  impulse/defect receipt is present. The articulated tree supports only fixed
  and scalar revolute/prismatic/helical joints.
- The articulated lane provides a full 6x6 spatial inertia and an articulated
  tree, but no geometry-derived inertia assembly, flexible body, closed-loop
  graph, or Euler-disc-specific rule. The separate point/impulse API remains
  parameterized by the legacy centre-of-mass diagonal principal inertia.
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
  momentum invariants hold to solve tolerance and roundoff for diagonal
  principal inertia; prescribed torque obeys the stated midpoint work balance.
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
patches to original COM/principal-frame `RigidBodyState` and `MassProperties`
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

During support, only the normal applied force is balanced. With tangent inverse mass `kappa = 1/m + radius^2/I`, rolling requires traction `-(F_t/m)/kappa` inside the circular friction limit. Collinear sliding retains analytic constant-force motion, exact slip-arrest times and resolved reversal. Frictionless support admits arbitrary tangent drive.

Resolved changing-direction sliding uses the original isotropic Coulomb disk in a velocity-implicit discretization, bounded by `max_sliding_step_s` and `max_sliding_substeps`. Actual translation, world momentum, contact-point slip and torque include the original spin inertia. Each accepted step uses its resolved constant force in the native trajectory; endpoint Coulomb feasibility does not assert an exact continuous force direction. The curved trajectory is first order and arrest is resolved at a timestep boundary. Refinement is necessary. Coarse trials that would make friction assist the initial slip are halved, including apparent endpoint arrest. Supercritical transverse drive cannot be silently converted to persistent rolling. Rejected trials count toward the work cap.

Both routes evaluate lost work from actual midpoint relative slip on each constant-force segment, never from the fictitious impulse-law loss at the force-driven free velocity. Curved-step work belongs to the declared time discretization and converges with its sliding timestep. Rolling traction has zero relative loss; the moving plane supplies signed work. No mechanical loss is automatically assigned to heat.

Attitude uses an exact exponential for commuting spin/torque and bounded second-order Lie midpoint substeps otherwise. Newly supplied penetration refuses. Each curved supported segment checks its complete quadratic normal-gap path against the original floating geometry allowance, in addition to actual endpoint kinematics. Private continuation retains only bounded contact roundoff; no pose or momentum is projected.

Every accepted advance independently checks total linear impulse, orbital-plus-spin angular impulse about the world origin, and `K_after - K_before = external_force_work + plane_work - dissipated_energy`. Cancellation, unresolved geometry/work, nonadvancing clocks and exhausted event, attitude or sliding budgets publish no partial state. Accumulating rebounds do not alter restitution.

Fifteen focused native integration groups pass: the eleven original ballistic, support, collinear, incline, spin, frame, attitude and refusal groups plus four curved-flow groups. An independent separable continuum solution for `u' = a - B*u/|u|` verifies velocity, position, spin and loss convergence at 16/32/64 steps, with error ratios 1.8–2.2; a separate 32/64/128 ladder verifies arrest followed by driven rolling. Coarse-step passivity, supercritical non-arrest, moving-plane work, split continuation, cancellation and late sliding-budget rollback are also exercised. These are focused actual-source native tests, not a full Cargo/DSR workspace gate. General multicontact Moreau-Jean dynamics, external torque, anisotropic bodies, multiple planes, finite-mass plane recoil, rolling/torsional resistance, arbitrary-shape certified CCD, source calibration and contact adjoints remain outside this owner.

## Native rigid-body simultaneous frictional impacts

`frictional_impact::apply_rigid_frictional_impacts` binds the simultaneous
normal/Coulomb owner to original rigid-body states. Actual masses, principal
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

A velocity-only event may contain dynamic bodies and prescribed supports sharing one common translation. Their total support work is that translation dotted with the uniquely determined net support impulse. Independently moving prescribed supports require a separate work-determinacy witness; this consumer returns `UnresolvedSupportWork` instead of publishing an arbitrary work allocation. Existing actual endpoint-velocity, kinetic/work, linear/angular-momentum and cancellation gates still precede atomic body publication.

Three new native G1/G3 groups in `tests/frictional_impact/redundant.rs` cover a finite-height, three-point body in stick and gross slide with independently derived load transfer, energy and momenta; frictionless duplicate normal rows with restitution; common conveyor translation/work and refusal of independently moving supports. Canonical input permutation is exact. These and four original impact regressions pass, together with four existing configuration-step regressions, in the bounded native source-selected lane. Generic rank/feedback, configuration, contact discovery and physical heat-partition boundaries remain explicit.


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
Children use the declared shared principal-axis convention; incompatible
geometry or nonrepresentable bookkeeping is refused.

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
