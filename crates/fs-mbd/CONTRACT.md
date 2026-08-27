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
a 15-DoF lower-body-and-waist reduction of Unitree's legacy `g1_29dof` URDF,
and the seven-axis KUKA LBR iiwa 7 R800 description from `iiwa_stack`. It adds
no second pose, twist, wrench, inertia, joint, or articulated-model type.

Canonical poses, twists, wrenches, adjoints, and coadjoints remain owned by
`fs-ga`; the articulated lane consumes those types rather than creating a
parallel robot-math representation. Contact, impacts, loop constraints,
free-floating-base equilibrium, and time integration remain dedicated owner
boundaries and are not approximated here.

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
  scalar-DoF indexing. The frame above the root is prescribed by `BaseState`,
  including its body twist and non-gravity body acceleration; it is not a
  solved floating body.
- `articulated::forward_kinematics` returns world-from-link `Se3` poses and
  body-coordinate twists. `inverse_dynamics` implements a recursive
  Newton-Euler pass and reports required generalized effort. `forward_dynamics`
  implements Featherstone's articulated-body algorithm without constructing a
  dense generalized mass matrix; supplied actuator efforts are checked against
  declared limits.
- `robot_models::unitree_g1_lower_body_waist_15dof` builds a fixed-pelvis tree
  with the source joint order of six left-leg, six right-leg, and three waist
  DoFs. It uses the link inertias, joint origins, axes, and hard limits from the
  official Unitree `g1_29dof.urdf` pinned in the returned metadata. The pinned
  Unitree README marks that legacy source variant deprecated.
- `robot_models::kuka_lbr_iiwa7_r800` builds the fixed-base seven-axis chain
  from `iiwa_stack`'s pinned `iiwa7.xacro`. Its 300 N m effort and 10 rad/s
  velocity limits are the Xacro's generic macro defaults, not claimed KUKA
  hardware limits.
- `robot_models::CatalogRobotModel` binds a validated `ArticulatedModel` to its
  stable compact joint order and immutable `RobotModelMetadata`. The metadata
  records exact pinned URLs, revisions, Git blob identities, source status,
  units, derivation, and material omissions. `ROBOT_MODEL_CATALOG_VERSION`
  versions this typed in-source layout; it is not a serialized URDF schema.

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
- Catalog builders deterministically preserve their declared root-first link
  topology and compact source-joint order. Their retained numeric records use
  URDF SI units and pass the same articulated topology, limit, pose, and
  physical-inertia validation as caller-built models.

## Error model

`DynamicsError` is the refusal channel for non-finite input or derived result,
unrepresentable finite magnitude, invalid mass, invalid or physically
inconsistent principal inertia, invalid orientation, invalid directional
effective mass, unsupported non-centre-of-mass reference point, and invalid
duration. There are no panics in production paths for those invalid inputs.
Arithmetic overflow or a non-finite attitude update, kinematic transform,
effective-mass query, or impulse-work result is returned through the checked
boundary; already completed earlier steps of `advance` remain committed.

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

The `robot_models` module additionally checks both catalogs' link/DoF counts
and stable source order, bilateral G1 neutral-origin symmetry, an independent
iiwa neutral endpoint, admitted hard limits and physical inertias, zero dense
generalized-matrix entries, zero-input ABA, deterministic rebuilds, and retained
provenance/omission records.

These are local G0-style checks. They do not constitute contact, constrained or
floating-base validation, full robot-model validation, performance evidence,
G4 fault injection, or G5 cross-ISA evidence.

## No-claim boundaries

- No contacts, impacts, friction, complementarity, penetration handling, or
  static-friction capacity are represented.
- No holonomic loop constraints, nonholonomic rolling constraints, RATTLE
  projection, generalized-alpha lane, Newton/Krylov solve, or constraint
  impulse/defect receipt is present. The articulated tree supports only fixed
  and scalar revolute/prismatic/helical joints.
- The articulated lane provides a full 6x6 spatial inertia and an articulated
  tree, but no geometry-derived inertia assembly, flexible body, closed-loop
  graph, free-floating-base solve, or Euler-disc-specific rule. The separate
  point/impulse API remains parameterized by the legacy centre-of-mass diagonal
  principal inertia.
- Catalog entries are transcriptions and explicit reductions, not complete or
  manufacturer-validated digital twins. The G1 entry omits all arm DoFs and
  links, fixed body/sensor attachments, contact geometry, and their mass rather
  than lumping omitted mass into retained links. The iiwa entry omits the world
  link, fixed massless flange, damping, soft safety limits, and payload/tool.
  Both use prescribed fixed-base semantics: neither solves a floating base.
- No catalog entry loads meshes, supplies collision/contact geometry, models
  actuators/gearing/transmissions, certifies hardware safety envelopes, parses
  URDF/Xacro at runtime, or claims the omitted upstream structures have no
  physical effect.
- No collision detection, signed gap, support mapping, common-point proof,
  contact selection, impact/restitution law, complementarity solve, friction
  cone, or no-slip constraint is implemented. An equal-and-opposite impulse is
  algebraic action/reaction only; it does not by itself establish angular-
  momentum conservation about a shared contact point or physical admissibility.
- The midpoint update is not claimed symplectic, variational, energy exactly
  conserving, momentum exactly conserving for arbitrary inertia, adaptive,
  adjoint-capable, or physically validated.
- Diagnostics describe the simulated smooth state; they do not mint a
  certificate, authority, or release-level conservation claim.
