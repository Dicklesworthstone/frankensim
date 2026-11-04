# Source-bound nonlinear rod motion

The opt-in `fs-couple/rod-dynamics` feature connects the source-bound Cosserat
rod, existing full-inertia validation and principal-frame conversion, and the
original `fs-time::rigid_body_momentum_step`. No material law, geometrical strain
map, analytical force calculation, eigensolver or free-rotation integrator is
reimplemented.

```sh
cargo test -p fs-solid --lib rod::material_field::motion
cargo test -p fs-solid --test rod_motion_model
cargo test -p fs-couple --features rod-dynamics --lib rod_motion::inertia::tests
cargo test -p fs-couple --features rod-dynamics --test rod_motion
cargo run -p fs-couple --features rod-dynamics --example rod_motion
cargo run -p fs-couple --features rod-dynamics --example rod_motion -- --eccentric --refinement 1
cargo run -p fs-couple --features rod-dynamics --example rod_motion -- --eccentric --refinement 2
cargo run -p fs-couple --features rod-dynamics --example rod_motion -- --eccentric --refinement 4
```

The example uses manufactured, not calibrated, material coefficients. Its
`--eccentric` mode adds a pre-tensioned link at rotating root/tip offsets.
`--refinement n` accepts integers in 1..=32 and reduces the timestep while
retaining the same 0.02-second interval. Standard output is CSV containing
energies, signed energy change, step balance residuals and material-tip
coordinates. Mode, timestep and inertia-lowering residual go to standard error.
These are reproduction commands, not a claim that native validation has run.

## One mechanical owner, two explicit inertia views

`SpatialMaterialRod::into_motion_model` retains the original centered-inertia
API. `into_spatial_motion_model` retains eccentric attachment inventory. Both
specialize the same `PreparedRodMotion<I>` and share source admission, frozen
thermal states, intrinsic geometry, stationary supports and native analytical
force evaluation. Root clamping or free release is explicit. A `RodMotionPose`
always contains MATERIAL-node world positions and canonical body-to-world
`So3` frames; it is never silently replaced by mass-centroid geometry.

Reference density and original annular geometry supply endpoint row-lumped
mass and section rotary inertia: `L0/2 sum(rho A)` and
`L0/2 sum(rho I) * diag(2,1,1)` at each segment endpoint. Each tensile link adds
half its original reference mass at each attachment. Stretching or heating does
not re-evaluate that inventory on the current length or density.

The centered constructor still rejects nonzero attachment offsets. The spatial
view retains actual nodal mass `M`, material-frame centroid `c`, and the complete
central tensor `I_C`, including products of inertia. Successive central-inertia
merges avoid subtracting large parallel-axis tensors. Symmetric eccentric
attachments can have zero centroid but still add substantial rotary inertia;
zero centroid does not make that mass equivalent to centered point mass.

This is the declared ENDPOINT lumping of strand mass, not distributed cable
kinetic energy, sag, internal cable modes or three-dimensional beam inertia.

## Shared dynamic trajectory

`RodMotion::new` initializes the centered view. `SpatialRodMotion::new_spatial`
initializes the full-inertia view and additionally requires a principal-tensor
reconstruction tolerance in `(0, 1e-8]`. The shared `SpatialInertia` validator
checks physical tensor consistency. Non-diagonal tensors use the existing
`fs-mbd-material::principal::PrincipalFrame` lowering once at initialization;
already-diagonal tensors retain their material axes without an eigensolve.
`maximum_inertia_tensor_residual()` exposes the measured lowering residual,
not an experimental uncertainty or an error certificate for the simulation.

Both views use the same kick/drift/kick and accounting implementation. State
stores total world linear momentum and CENTRAL spin in the original material
axes. Input and recovered velocities refer to the material node. With material
orientation `R` and body angular velocity `omega`, initialization uses

```text
x_C = x_O + R c
v_C = v_O + R (omega cross c)
p   = M v_C
L_C = I_C omega
```

The kinetic drift advances the centroid and the native gyroscopic rotor, then
updates material position without subtracting large absolute centroid positions:

```text
x_O_next = x_O + dt p/M + R c - R_next c
```

Principal axes are only a fixed computational basis; they never change material
directors, source orientation or reference geometry. Cloning shares immutable
sources and kinetic frames. Each accepted step reuses its endpoint force for
the next start, requiring one new analytical force assembly, no dense tangent
and no eigendecomposition in the step.

## Wrenches, work and support reactions

`RodNodalWrench` means a force AT the material node plus a pure couple, both
held in world coordinates for the interval. Native analytical torques are also
about the material node. The dynamic owner shifts them to the actual centroid:

```text
tau_C_body = tau_O_body - c cross (R^T F)
```

Applied work uses centroid force power plus that shifted couple power. This is
identically the material-node force/couple power; omitting the shift would
create false work under eccentric rotation. Energy includes kinetic storage
and native elastic/link/gravity potential. Gravity is not added again as source
work, but its offset moment participates in external angular-impulse accounting.
Total angular momentum is centroid orbital momentum plus central spin.

Free eccentric nodes and fixed material frames are supported. A fixed frame
may translate in its free world components; its reported pure support couple
ABOUT THE MATERIAL NODE includes the eccentric translational reaction. Fully
clamped nodes remain stationary. Centered nodes retain their original component
position supports. An eccentric node with free rotation and ANY constrained
material-position component returns `EccentricPositionSupport`: its pin/slider
requires coupled constraint dynamics, not independent momentum projection.

## Limits and verification status

This remains second-order fixed-step symplectic splitting on smooth conservative
paths, **not** an exactly energy-conserving discrete-gradient method. Numerical
energy defects are reported and bounded before candidate publication; they are
neither physical heat nor exact continuous-trajectory work. An energy gate is
not a stability or accuracy proof. Elastic-wave timestep and mesh refinement
remain necessary. Swept-angle limits do not certify self-contact avoidance or
rotation-logarithm branch avoidance between samples.

Material frequency restrictions, fixed thermal/source states, 64-segment and
link budgets remain enforced. Explicit angular material damping is described
below; plasticity, contact, moving supports, live material replacement and
distributed cable inertia remain outside this path. Failed
or cancelled steps preserve the complete accepted state; clones retain in-process
restart state. Principal-axis selection is not differentiated through repeated
eigenvalues.

The eccentric extension adds six source-inventory and eight trajectory Rust
regression groups, including direct constituent energy, a full 6x6 spatial-mass
oracle, eccentric frame reactions, free fall, free-flight refinement, centered
specialization parity, and rejection/restart. At implementation time these
native tests were not executed because the available environment lacked
rustc/Cargo. Independent Python equation checks covered 1,000 inventory tensors,
4,000 constituent kinetic energies and 4,000 spatial kick/reaction/drift cases.
An independently differenced nonlinear rod showed fourfold endpoint energy-drift
reduction per timestep halving. Those checks are not native Rust build proof,
spatial-convergence evidence or calibrated material validation.

## Source-bound bending and torsional damping

`PreparedRodBendingLoss` connects each original circular section's explicit
`ConcentricKelvinVoigtSection` to the native trajectory. The material owner
integrates the layer viscosities as `C_b = sum(eta I)`, with no averaged material
or fitted damping ratio. Its two bending resistances are `C_b/L0`.
`PreparedRodTorsionalLoss` independently supplies `sum(eta_s J)/L0`; bending
viscosity does not authorize torsional viscosity. Both source bundles must retain
the actual rod's geometry, original material cards, full physical query and
mechanical-property selections. Transient loss requires explicit broadband
support at the held state. After a material/temperature refresh, prepare new
matching loss sections before stepping.

```rust,ignore
use fs_couple::rod_motion::{PreparedRodBendingLoss, RodAngularLosses};

let bending = PreparedRodBendingLoss::new(state.model(), bending_sections, &gate)?;
let losses = RodAngularLosses { bending: Some(&bending), torsional: None };
let (candidate, receipt) = state.step_damped_candidate(
    dt, losses, &world_wrenches, controls, &gate,
)?;
let physical_bending_heat = receipt.bending_heat.as_ref().unwrap().total_heat_j;
```

`step_spatial_damped_candidate` accepts the existing `PreparedRodLoads` for
rotating material-point forces. Either inertia view admits the same loss API.
The existing torsional-only step APIs retain their original flow.

For `phi = Log(R_i^-1 R_j)`, bending uses `d(phi_y)/dt` and `d(phi_z)/dt`.
Exact SO(3) Jacobian pullbacks make those rates work-conjugate to the nodal
moments under finite rotations; a common rigid spin produces no material loss.
Each frozen-pose scalar flow relaxes exponentially through the original central
inertia and reports its positive integrated work. Supported frames receive
the corresponding reaction. Palindromic sweeps and the surrounding symmetric
loss/elastic/loss split retain second-order temporal convergence. Strong
viscosity adds no explicit damping stability restriction; elastic wave resolution
and timestep refinement still matter.

`RodMotionStep::bending_heat` and `torsional_heat` report separate total and
original-segment physical heat. The signed integration defect includes both
exactly once. This fixed-state mechanical path reports heat without depositing
it into a thermal reservoir. This angular path supplies no axial/shear strain viscosity,
hereditary bending spectrum or thermoelastic bending feedback. An omitted law
is unmodeled, and an explicit zero-viscosity law leaves the native trajectory
unchanged. Focused regressions exercise independent damped-oscillator and heat
refinement, both bending planes, stiff viscosity, fixed supports, finite-rotation
objectivity with eccentric full inertia, and atomic source/refusal behavior.

## Source-bound axial damping

`PreparedRodAxialLoss` connects the existing `ConcentricAxialViscosity` material
owner to the same trajectory. Each original segment retains its actual layered
section, material receipts, state and geometry. The resistance is
`C_a/L0 = sum(eta_a A)/L0`; the viscosity must explicitly support a broadband
Kelvin law at the held source state. A frequency-local loss sample or a missing
viscosity property cannot supply that law.

```rust,ignore
use fs_couple::rod_motion::{PreparedRodAxialLoss, RodAngularLosses};

let axial = PreparedRodAxialLoss::new(state.model(), axial_sections, &gate)?;
let angular = RodAngularLosses { bending: None, torsional: None };
let (candidate, receipt) = state.step_axially_damped_candidate(
    dt, &axial, angular, &world_wrenches, controls, &gate,
)?;
let physical_axial_heat = receipt.axial_heat.as_ref().unwrap().total_heat_j;
```

`step_spatial_axially_damped_candidate` accepts prepared material-point loads.
Either entry point can include independently sourced bending and torsion.
The rate is the derivative of the native Cosserat extension `L0*gamma_x`,
including its geodesic-midpoint frame dependence. The force and torque are
therefore conjugate to the same finite-rotation strain used by elastic motion;
common rigid translation and spin do not become axial viscosity. The pullback
retains eccentric centers of mass, full central inertia and stationary support
reactions. Each frozen-pose scalar relaxation uses an exact exponential impulse
and nonnegative physical heat. Bounded stack scratch and symmetric segment and
mechanism sweeps avoid a dense damping solve and retain second-order splitting.
Elastic-wave resolution and timestep refinement remain necessary.

`RodMotionStep::axial_heat` reports total and original-segment heat independently
of bending, torsion and the signed integration defect. The defect includes each
physical loss once. This fixed-state entry point exports heat without depositing
it into a thermal store, and supplies no transverse-shear viscosity or hereditary
axial memory. Explicit zero axial viscosity preserves the existing trajectory.
Source mismatch, cancellation and a failed energy gate leave the accepted state
unchanged. The focused integration fixtures cover axial ringdown and heat
refinement, stiff viscosity, reactions, finite-rotation objectivity, zero-law
parity and source/refusal behavior.
