# CONTRACT: fs-versatran-wasm

## Purpose and layer

Layer: **L6 browser boundary**. This standalone crate exposes the
source-bounded AMF Versatran composition for US 3,212,649. The generic joint
owner is `fs_mbd::versatran`; this package only admits browser scalar inputs
and serializes a stable envelope.

Claim 1 identifies six actuator-controlled motions: column rotation; vertical
and horizontal movement of the horizontal arm; wrist rotation about a
horizontal axis; wrist swing about a central vertical axis; and operation of
a work-manipulating member. The first five are member-motion joints. The sixth
is represented only as the disclosed reciprocating work-member rack/sleeve
operating coordinate, not as a sixth rigid-pose degree of freedom.

## Public surface

`versatran_topology_step(column_rotation_rad, arm_vertical_normalized,
arm_horizontal_normalized, wrist_rotation_rad, wrist_swing_rad,
work_member_rack_normalized, automatic_program_mode_selected) -> String`
returns either:

- `ok`: a six-scalar generic-joint receipt with five geometric member-motion
  joints, one internal work-member operation channel, three revolute axes,
  three prismatic axes, the admitted coordinates and automatic-mode predicate;
  or
- `refusal`: a stable code, diagnostic, and repair.

## Invariants

1. The composition is exactly three revolute and three prismatic scalar joints.
2. Only five coordinates move the named column, arm, or wrist members. The
   sixth operates the work-manipulating member internally and cannot be
   described as a six-axis rigid pose.
3. The vertical and horizontal arm coordinates plus work-member rack coordinate
   are dimensionless normalized presentation values in `[0,1]`, not metres.
4. `+Y` denotes the claim's vertical direction and `+X` the horizontal
   arm/rack direction solely in the museum normalization convention. They are
   not original AMF dimensions, zero transforms, or a claimed global frame.
5. Selecting automatic-program mode reports a claimed control-mode selection;
   it does not run a program, interpolate a trajectory, or manufacture
   feedback/controller behavior.
6. The result explicitly states that historical geometry and dynamics are not
   available from this source boundary.

## Error model

`input-outside-domain` refuses non-finite rotations or normalized arm/rack
coordinates outside `[0,1]`. `multibody-refusal` preserves an unexpected
generic-joint-owner refusal. Every refusal includes a concise explanation and
repair.

## Determinism class

This is a pure synchronous query. Equivalent admitted inputs produce
byte-identical output on the same Rust/WASM build. It has no clock, randomness,
host identity, stored program, or hidden mutable state.

## Cancellation behavior

The bounded one-state query starts no task and performs no I/O, so it has no
cancellation seam.

## Unsafe boundary

None. This crate forbids unsafe Rust.

## Feature flags

None. `wasm-bindgen` is target-selected for `wasm32`; native tests exercise the
same pure boundary function.

## Conformance

Native tests cover the three-revolute/three-prismatic composition, the
five-geometric-plus-one-work-operation distinction, automatic-mode reporting,
and typed refusal. `fs-mbd` separately checks generic-joint construction and
input admission.

## No-claim boundary

This package does not provide machine dimensions, transforms, forward or
inverse kinematics, collision geometry, link mass/inertia, cylinder stroke,
pressure, flow, force, torque, speed, payload, cycle timing, accuracy,
controller gains, a recorded program, actuator feedback, or commercial/safety
performance. It neither reconstructs gripper jaws nor claims calibrated AMF
hardware behavior.
