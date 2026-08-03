# CONTRACT: fs-mbd

## Purpose and layer

Layer: **L3 FLUX**. `fs-mbd` currently provides the safe-Rust foundation for
one smooth, unconstrained rigid body whose reference point is its centre of
mass. It is intentionally a bounded dependency-free core: validated diagonal
principal inertia, body-to-world orientation, momentum-form state, uniform
gravity, a constant external wrench, diagnostics, and a deterministic fixed
step update.

This is a foundation for later rigid multibody work. It does not substitute for
the blocked contact, joint, constraint, nonlinear-solver, or nonholonomic
lanes.

## Public types and semantics

- `Vec3` is a small Cartesian vector. Every public field that carries a vector
  names its frame in the field name or documentation.
- `UnitQuaternion` is a normalized body-to-world quaternion. Construction
  rejects zero/non-finite norm, normalizes, and chooses the representative for
  which the first nonzero component of `(w, x, y, z)` is positive. This makes
  the quaternion double cover deterministic.
- `Pose` combines a world-space centre-of-mass position with that orientation.
- `MassProperties::new(mass, center_of_mass_body, principal_inertia_body)`
  accepts finite positive mass and principal moments satisfying the rigid-body
  triangle inequalities. The centre of mass must be exactly zero in the
  principal frame; offset-reference spatial inertia is refused.
- `RigidBodyState` stores world-frame linear momentum and body-frame angular
  momentum. It has no independent velocity or angular-velocity state.
- `Wrench` carries a constant world-frame force and principal-body-frame torque
  for one step. `Gravity` is a uniform world-frame acceleration.
- `RigidBodyIntegrator::step` returns a `StepReceipt` with the before/after
  states and diagnostics. Translation uses the midpoint linear momentum under
  the constant total force. Body angular momentum uses midpoint RK2 for
  `Ldot = L × I⁻¹L + tau`; attitude uses the body midpoint angular velocity in
  a right-composed quaternion exponential update.
- `DynamicsDiagnostics` reports translational/rotational kinetic energy,
  gravity potential with zero at the world origin, total mechanical energy,
  world linear momentum, and world angular momentum about that origin.
- `RigidBodyIntegrator::advance` invokes its cancellation callback before each
  whole step. It returns `AdvanceOutcome::Cancelled` with the last fully
  committed state and diagnostics when cancellation is observed.

## Invariants

- Public construction rejects non-finite mass-property, state, gravity, wrench,
  duration, and quaternion inputs at the relevant boundary.
- The public orientation remains normalized and uses a canonical quaternion
  sign after every successful attitude update.
- The declared force and gravity are world-frame quantities; the declared
  torque and angular momentum are principal-body-frame quantities. Diagnostics
  rotate the latter to world space before adding orbital angular momentum.
- Each successful `step` constructs and returns a complete new state. In
  `advance`, only complete preceding steps are committed when cancellation is
  observed.

## Error model

`DynamicsError` is the refusal channel for non-finite input, invalid mass,
invalid or physically inconsistent principal inertia, invalid orientation,
unsupported non-centre-of-mass reference point, and invalid duration. There
are no panics in production paths for those invalid inputs. Arithmetic overflow
or a non-finite attitude update is returned through the same checked state and
orientation construction; already completed earlier steps of `advance` remain
committed.

## Determinism class

The step order, frame transformations, quaternion sign convention, and
cancellation boundary are deterministic functions of the admitted inputs. The
implementation currently uses the Rust standard library's floating-point
`sqrt`, `sin`, and `cos`; therefore this contract makes no cross-ISA bitwise
claim and retains no cross-ISA evidence. Same-platform repeatability is an
implementation property to be exercised by consumers, not a retained G5
promotion from this seven-test foundation.

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

None. This initial core depends only on `std`/`core` and has no feature flags.

## Focused evidence

The inline test module in `src/lib.rs` contains seven focused checks:

1. mass, inertia, and centre-of-mass-reference refusals;
2. canonical quaternion double-cover selection and a known rotation;
3. analytic constant-gravity translation plus energy conservation;
4. constant force/torque updates in their declared frames;
5. torque-free spherical-body energy and world angular-momentum conservation;
6. refinement reducing the measured energy defect for an asymmetric free body;
7. cancellation before a whole step, preserving the last completed state.

These are local G0-style checks. They do not constitute full multibody,
constraint, contact, physical-validation, performance, G4 fault-injection, or
G5 cross-ISA evidence.

## No-claim boundaries

- No contacts, impacts, friction, complementarity, penetration handling, or
  static-friction capacity are represented.
- No joints, holonomic constraints, nonholonomic rolling constraints, RATTLE
  projection, generalized-alpha lane, Newton/Krylov solve, or constraint
  impulse/defect receipt is present.
- No full 6x6 spatial inertia, non-principal inertia tensor, offset reference
  point, material receipt, geometry-derived inertia assembly, flexible body,
  articulated multibody graph, or Euler-disc-specific rule is present.
- The midpoint update is not claimed symplectic, variational, energy exactly
  conserving, momentum exactly conserving for arbitrary inertia, adaptive,
  adjoint-capable, or physically validated.
- Diagnostics describe the simulated smooth state; they do not mint a
  certificate, authority, or release-level conservation claim.
