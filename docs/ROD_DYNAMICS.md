# Source-bound nonlinear rod motion

The opt-in `fs-couple/rod-dynamics` feature connects the existing source-bound
Cosserat rod to the original `fs-time::rigid_body_momentum_step`. No material
law, geometrical strain map, analytical force calculation or free-rotation
integrator is duplicated.

```sh
cargo test -p fs-solid --test rod_motion_model
cargo test -p fs-couple --features rod-dynamics --test rod_motion
cargo run -p fs-couple --features rod-dynamics --example rod_motion
```

The example uses explicitly manufactured coefficients, writes CSV physical
observables, and is not a calibrated material or instrument demonstration.

## Mechanical model

`SpatialMaterialRod::into_motion_model` freezes the original source bundles,
reference geometry, thermal state and supports. Root clamping or free release
is explicit. A `RodMotionPose` contains the existing world positions and
canonical body-to-world `So3` frames. The native analytical owner evaluates
full nodal forces and body couples, including root reactions, without allocating
a dense tangent. Source bounds and rotation-logarithm branch cuts still apply.

Reference density and original annular geometry supply endpoint row-lumped
mass and cross-section rotary inertia: `L0/2 sum(rho A)` and
`L0/2 sum(rho I) * diag(2,1,1)` at each segment endpoint. Orbital angular momentum
is computed separately as `(x-origin) cross p`. Heating or stretching never
re-evaluates inventory on current geometry. Centerline tensile links add half
their reference mass at each endpoint. Eccentric attachments refuse because
this diagonal model cannot represent their coupled spatial inertia.

`RodMotion` stores world linear and body angular momenta. A step composes half
force/couple kicks, translational drift, native gyroscopic rotation, and the
second kick. Its accepted endpoint force is reused as the next starting force;
only one new analytical force assembly is needed per step. Immutable source
bundles are shared when the trajectory is cloned. Stationary supported degrees
of freedom do not drift, and their balancing impulses are reported.

## Accounting and limits

Energy includes kinetic storage and native elastic/link/gravity potential.
Applied world-frame wrenches use symmetric kick power quadrature; gravity is
not counted twice as applied work. The signed energy defect is reported and
bounded before publishing a candidate. It is neither physical loss nor heat,
nor an exact continuous-trajectory work or error certificate. Momentum reports
include applied loads, gravity, support reactions and orbital plus section-spin
terms about the declared origin.

This is a second-order fixed-step symplectic splitting on smooth conservative
paths, **not** the exactly energy-conserving discrete-gradient reference sought
by the full dynamic-rod bead. Explicit elastic-wave stability and timestep/mesh
refinement remain necessary. An energy-defect gate alone is not a stability or
accuracy proof. The swept-angle cap prevents excessive per-node rotation but
does not certify absence of self-contact or branch crossings between samples.
Finite source frequency bands refuse because arbitrary transient motion is not
proved band-limited by timestep selection. No material damping, plasticity,
contact, moving supports, live source refresh, variable cross-section stretch,
or distributed cable-inertia claim is introduced. Existing 64-segment and
source budgets remain unchanged. Refusal or cancellation leaves the complete
accepted motion reusable; cloning provides deterministic in-process restart.
