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
Contact, impacts, loop constraints, and time integration remain dedicated owner
boundaries and are not approximated here. Control policy is likewise a separate
owner: this crate consumes declared efforts and external wrenches but invents
neither.

The `goddard` module also exposes a narrow museum-facing composition for US
1,102,653. It reuses the single-body owner for torque-free primary and
gyroscope display poses, then reports only source sequence predicates and the
Claim 2 tapered-tube ratio. The older `step_goddard_rocket` liquid-nozzle
calculation remains an explicitly adjacent interpretive model; it is not a
model of the mechanism claimed in US 1,102,653.

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
  the constant total force. Body angular momentum uses midpoint RK2 for
  `Ldot = L × I⁻¹L + tau`; attitude uses the body midpoint angular velocity in
  a right-composed quaternion exponential update.
- `RigidBodyIntegrator::diagnostics` returns
  `Result<DynamicsDiagnostics, DynamicsError>`. It revalidates its state,
  gravity, and mass properties and refuses non-finite derived energies or world
  angular momentum rather than publishing an `inf`/`NaN` ledger.
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
duration. There are no panics in production paths for those invalid inputs.
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
9. refinement reducing the measured energy defect for an asymmetric free body;
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
  projection, generalized-alpha lane, Newton/Krylov solve, or constraint
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
- The midpoint update is not claimed symplectic, variational, energy exactly
  conserving, momentum exactly conserving for arbitrary inertia, adaptive,
  adjoint-capable, or physically validated.
- Diagnostics describe the simulated smooth state; they do not mint a
  certificate, authority, or release-level conservation claim.
- The source-bounded Goddard apparatus lane claims no absolute patent
  dimensions, mass properties, burn rate, force, thrust, Mach number,
  trajectory, aerodynamic stability, gyroscope torque capacity, or empirical
  validation. Its ideal zero-world-rate camera support is a kinematic teaching
  limit, not a prediction that the printed mechanism rejects arbitrary loads.
- `goddard::step_goddard_rocket` uses an adjacent liquid-propellant/de Laval
  interpretation. No consumer may present those outputs as a claim, disclosed
  embodiment, or source-derived parameter of US 1,102,653.
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
